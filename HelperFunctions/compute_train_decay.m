function train_results = compute_train_decay(VibData, do_plot)
% COMPUTE_TRAIN_DECAY  Within-train dynamics of vibrotactile burst trains.
%
%   train_results = compute_train_decay(VibData)
%   train_results = compute_train_decay(VibData, do_plot)
%
%   For each train, extracts the Hilbert envelope over the full train window
%   and, per burst, measures 10–90% rise time, the inter-burst decay time
%   constant (log-linear fit over the 85%->40% band), the residual amplitude
%   at the next burst onset (% of peak), and — for the final burst — the
%   settling time to 5% of peak with an adaptive noise floor.
%
%   Output struct array fields (one row per burst within each train):
%     .Position .PotStep .IBI_ms .RepIdx .BurstSeq — condition identifiers
%     .RiseTime_ms     — 10–90% rise time (ms)
%     .Tau_ms          — decay time constant after burst offset (ms)
%     .Settle_Meas_ms  — settling time, final burst (ms)
%     .Residual_Pct    — envelope at next onset as % of peak (0 for final burst)
%     .Peak_g          — envelope peak of the burst (g)

    %% === PARAMETERS ===
    if nargin < 2, do_plot = true; end
    FIT_UPPER = 0.85; 
    FIT_LOWER = 0.40;
    SETTLE_TARGET = 0.05; % 5% of peak for final settle
    
    %% === INITIALIZATION ===
    train_results = struct('Position',{}, 'PotStep',{}, 'IBI_ms',{}, ...
                           'RepIdx',{}, 'BurstSeq',{}, 'RiseTime_ms',{}, 'Tau_ms',{}, ...
                           'Settle_Meas_ms',{}, 'Residual_Pct',{}, 'Peak_g',{});
    idx = 0;
    
    % For diagnostic plotting of whole trains
    examples = struct();
    ex_count = 0;
    seen_combos = {}; 

    fprintf('\n================================================================\n');
    fprintf('  TRAIN DYNAMICS (Rise, IBI Decay, Final Settle)\n');
    fprintf('  v3 protocol: 3 pot x 4 IBI x 5 reps, 6 bursts/train\n');
    fprintf('================================================================\n');

    nPos = numel(VibData);
    
    for p_idx = 1:nPos
        D = VibData(p_idx).Raw;
        if isempty(D), continue; end
        
        Fs = D.Fs;
        dt = 1/Fs;
        
        % 20Hz HPF to remove gravity/drift
        [b, a] = butter(4, 20/(Fs/2), 'high');
        Z_filt = filtfilt(b, a, double(D.Z));
        
        if ~isfield(VibData(p_idx), 'Trains') || isempty(VibData(p_idx).Trains)
            continue;
        end
        Trains = VibData(p_idx).Trains;
        
        for t_i = 1:numel(Trains)
            T = Trains(t_i);
            if T.NBursts == 0, continue; end
            
            % --- 1. EXTRACT FULL TRAIN WINDOW ---
            burst_len_samps = round((T.BurstMs / 1000) * Fs);
            first_onset = T.BurstOnsetIdx(1);
            last_onset  = T.BurstOnsetIdx(end);
            last_fall   = last_onset + burst_len_samps;
            
            pre_pad  = round(0.1 * Fs); % 100ms before train
            post_pad = round(0.8 * Fs); % 800ms after train
            
            w_start = first_onset - pre_pad;
            w_end   = last_fall + post_pad;
            if w_start < 1 || w_end > numel(Z_filt), continue; end
            
            raw_train = Z_filt(w_start : w_end);
            try env_raw = abs(hilbert(raw_train)); catch, env_raw = abs(raw_train); end
            
            % Smooth envelope slightly (4ms window)
            sm_win = round(0.004 * Fs);
            env_train = movmean(env_raw, sm_win);
            
            % Time axis for the whole train (0 = first onset)
            t_axis_train = ((1:numel(raw_train)) - pre_pad - 1) * dt * 1000;
            
            train_peaks  = zeros(1, T.NBursts);
            train_rises  = zeros(1, T.NBursts);
            train_taus   = NaN(1, T.NBursts);
            train_resids = NaN(1, T.NBursts);
            train_settle = NaN(1, T.NBursts);
            
            % --- 2. PROCESS EACH BURST IN THE TRAIN ---
            for b_seq = 1:T.NBursts
                b_onset_abs = T.BurstOnsetIdx(b_seq);
                b_fall_abs  = b_onset_abs + burst_len_samps;
                
                % Convert absolute indices to snippet-relative indices
                idx_onset = b_onset_abs - w_start + 1;
                idx_fall  = b_fall_abs - w_start + 1;
                
                % Clamp to valid range
                idx_fall = min(idx_fall, numel(env_train));
                
                % A. FIND LOCAL PEAK
                burst_win = env_train(idx_onset : idx_fall);
                [local_peak, peak_loc] = max(burst_win);
                idx_peak = idx_onset + peak_loc - 1;
                train_peaks(b_seq) = local_peak;
                
                if local_peak > 1e-4
                    % B. RISE TIME (10% to 90%)
                    thresh_10 = local_peak * 0.10;
                    thresh_90 = local_peak * 0.90;
                    
                    rise_seg = env_train(idx_onset:idx_peak);
                    i_90 = find(rise_seg >= thresh_90, 1, 'first');
                    i_10 = find(rise_seg >= thresh_10, 1, 'first');
                    
                    if ~isempty(i_90) && ~isempty(i_10) && (i_90 > i_10)
                        train_rises(b_seq) = (i_90 - i_10) * dt * 1000; 
                    else
                        train_rises(b_seq) = NaN;
                    end
                    
                    % C. DECAY ANALYSIS
                    if b_seq < T.NBursts
                        % INTRA-TRAIN (Bursts 1 to N-1)
                        next_onset_abs = T.BurstOnsetIdx(b_seq+1);
                        idx_next_onset = next_onset_abs - w_start + 1;
                        
                        end_val = env_train(idx_next_onset);
                        train_resids(b_seq) = (end_val / local_peak) * 100;
                        
                        fit_seg = env_train(idx_fall : idx_next_onset);
                        i_upper = find(fit_seg < local_peak * FIT_UPPER, 1, 'first');
                        i_lower = find(fit_seg < local_peak * FIT_LOWER, 1, 'first');
                        
                        if ~isempty(i_upper) && ~isempty(i_lower) && (i_lower > i_upper)
                            y_fit = fit_seg(i_upper:i_lower);
                            t_fit = (0:(numel(y_fit)-1))' * dt;
                            valid = y_fit > 1e-9;
                            if sum(valid) > 3
                                pp = polyfit(t_fit(valid), log(y_fit(valid)), 1);
                                train_taus(b_seq) = -1 / pp(1) * 1000;
                            end
                        end
                        
                        % Settle inside train?
                        i_settle = find(fit_seg < local_peak * SETTLE_TARGET, 1, 'first');
                        if ~isempty(i_settle)
                            train_settle(b_seq) = i_settle * dt * 1000;
                        end
                        
                    else
                        % FINAL BURST
                        train_resids(b_seq) = 0; 
                        
                        fit_seg = env_train(idx_fall : end);
                        i_upper = find(fit_seg < local_peak * FIT_UPPER, 1, 'first');
                        i_lower = find(fit_seg < local_peak * FIT_LOWER, 1, 'first');
                        
                        if ~isempty(i_upper) && ~isempty(i_lower) && (i_lower > i_upper)
                            y_fit = fit_seg(i_upper:i_lower);
                            t_fit = (0:(numel(y_fit)-1))' * dt;
                            valid = y_fit > 1e-9;
                            if sum(valid) > 3
                                pp = polyfit(t_fit(valid), log(y_fit(valid)), 1);
                                train_taus(b_seq) = -1 / pp(1) * 1000;
                            end
                        end
                        
                        % Tail noise to prevent infinite wait on ringing
                        tail_start = max(1, numel(env_train) - round(0.1*Fs));
                        tail_noise = mean(env_train(tail_start:end)) + ...
                                     4*std(env_train(tail_start:end));
                        thresh = max(local_peak * SETTLE_TARGET, tail_noise);
                        
                        i_settle = find(fit_seg < thresh, 1, 'first');
                        if ~isempty(i_settle)
                            train_settle(b_seq) = i_settle * dt * 1000;
                        end
                    end
                end
                
                % Store per-burst
                idx = idx + 1;
                train_results(idx).Position       = VibData(p_idx).Name;
                train_results(idx).PotStep        = T.PotStep;
                train_results(idx).IBI_ms         = T.IBI_ms;
                train_results(idx).RepIdx         = T.RepIdx;
                train_results(idx).BurstSeq       = b_seq;
                train_results(idx).RiseTime_ms    = train_rises(b_seq);
                train_results(idx).Tau_ms         = train_taus(b_seq);
                train_results(idx).Settle_Meas_ms = train_settle(b_seq);
                train_results(idx).Residual_Pct   = train_resids(b_seq);
                train_results(idx).Peak_g         = local_peak;
            end
            
            % --- 3. SAVE FULL TRAIN EXAMPLE (rep 1 of each pot x IBI, first position) ---
            if T.RepIdx == 1
                combo_key = sprintf('%s_P%d_IBI%.0f', VibData(p_idx).Name, T.PotStep, T.IBI_ms);
                
                if ~ismember(combo_key, seen_combos)
                    seen_combos{end+1} = combo_key;
                    
                    ex_count = ex_count + 1;
                    examples(ex_count).Name    = VibData(p_idx).Name;
                    examples(ex_count).t_axis  = t_axis_train;
                    examples(ex_count).raw     = raw_train;
                    examples(ex_count).env     = env_train;
                    examples(ex_count).pot     = T.PotStep;
                    examples(ex_count).ibi     = T.IBI_ms;
                    
                    rel_onsets = (T.BurstOnsetIdx - w_start + 1) * dt * 1000 - pre_pad*dt*1000;
                    rel_falls  = rel_onsets + T.BurstMs;
                    
                    examples(ex_count).onsets = rel_onsets;
                    examples(ex_count).falls  = rel_falls;
                    
                    examples(ex_count).mean_rise  = mean(train_rises, 'omitnan');
                    examples(ex_count).mean_resid = mean(train_resids(1:end-1), 'omitnan');
                    examples(ex_count).final_tau  = train_taus(end);
                    examples(ex_count).final_set  = train_settle(end);
                    examples(ex_count).settle_array = train_settle;
                end
            end
        end
        fprintf('  %s: done\n', VibData(p_idx).Name);
    end
    
    %% === PLOT DIAGNOSTICS ===
    if do_plot && ex_count > 0
        plot_full_train_diagnostics(examples);
    end
    
    %% === SUMMARY TABLES ===
    if idx == 0, warning('No train results computed.'); return; end
    
    all_pot  = [train_results.PotStep];
    all_ibi  = [train_results.IBI_ms];
    all_seq  = [train_results.BurstSeq];
    all_rise = [train_results.RiseTime_ms];
    all_res  = [train_results.Residual_Pct];
    all_tau  = [train_results.Tau_ms];
    all_set  = [train_results.Settle_Meas_ms];
    
    u_pot = unique(all_pot);
    u_ibi = unique(all_ibi);
    n_bursts_per_train = max(all_seq);
    
    for pi = 1:numel(u_pot)
        pot = u_pot(pi);
        fprintf('\n----------------------------------------------------------------------------------\n');
        fprintf('  TRAIN DYNAMICS SUMMARY — Pot %d (250 Hz)\n', pot);
        fprintf('----------------------------------------------------------------------------------\n');
        fprintf('%-10s | %-12s | %-14s | %-14s | %-12s\n', ...
            'IBI (ms)', 'Rise (ms)', 'Resid %% (IBI)', 'Tau Final(ms)', 'Settle Final');
        fprintf('----------------------------------------------------------------------------------\n');
        
        for ii = 1:numel(u_ibi)
            ibi = u_ibi(ii);
            m_all = (all_pot == pot) & (all_ibi == ibi);
            m_ibi = m_all & (all_seq < n_bursts_per_train);  % intra-train bursts
            m_fin = m_all & (all_seq == n_bursts_per_train);  % final burst
            
            if ~any(m_all), continue; end
            
            fprintf('%8.2f   | %5.1f +/-%-4.1f | %5.1f%% +/-%-4.1f | %5.1f +/-%-5.1f | %5.0f +/-%-3.0f\n', ...
                ibi, ...
                mean(all_rise(m_all),'omitnan'), std(all_rise(m_all),'omitnan'), ...
                mean(all_res(m_ibi),'omitnan'),  std(all_res(m_ibi),'omitnan'), ...
                mean(all_tau(m_fin),'omitnan'),  std(all_tau(m_fin),'omitnan'), ...
                mean(all_set(m_fin),'omitnan'),  std(all_set(m_fin),'omitnan'));
        end
    end
    fprintf('----------------------------------------------------------------------------------\n');
end

%% === LOCAL PLOTTING: one figure per pot step, subplots by IBI ===
function plot_full_train_diagnostics(ex)
    all_pots = [ex.pot];
    unique_pots = sort(unique(all_pots));
    
    for pi = 1:numel(unique_pots)
        this_pot = unique_pots(pi);
        pot_mask = (all_pots == this_pot);
        pot_ex = ex(pot_mask);
        
        if isempty(pot_ex), continue; end
        
        % Sort by IBI ascending
        [~, sort_idx] = sort([pot_ex.ibi]);
        pot_ex = pot_ex(sort_idx);
        
        % Group by IBI — take first position only for cleaner plots
        seen_ibi = [];
        plot_ex = [];
        for ei = 1:numel(pot_ex)
            if ~ismember(pot_ex(ei).ibi, seen_ibi)
                seen_ibi(end+1) = pot_ex(ei).ibi; %#ok<AGROW>
                if isempty(plot_ex)
                    plot_ex = pot_ex(ei);
                else
                    plot_ex(end+1) = pot_ex(ei); %#ok<AGROW>
                end
            end
        end
        
        n_plots = numel(plot_ex);
        
        figure('Name', sprintf('Train Dynamics — Pot %d', this_pot), ...
               'Color', 'w', 'Position', [80+pi*30, 80+pi*30, 1400, 250*n_plots]);
        
        for i = 1:n_plots
            subplot(n_plots, 1, i); hold on; box on; grid on;
            E = plot_ex(i);
            
            % Plot signals
            h_raw = plot(E.t_axis, E.raw, 'Color', [0.85 0.85 0.85], ...
                         'LineWidth', 0.5, 'DisplayName', 'Raw Z');
            h_env = plot(E.t_axis, E.env, 'k-', 'LineWidth', 1.2, ...
                         'DisplayName', 'Envelope');
            
            h_on = []; h_off = []; h_set = [];
            n_bursts = numel(E.onsets);
            
            for j = 1:n_bursts
                % TTL ON / OFF lines
                if j == 1 && i == 1
                    h_on  = xline(E.onsets(j), 'g-', 'LineWidth', 1.5, 'DisplayName', 'TTL ON');
                    h_off = xline(E.falls(j),  'r-', 'LineWidth', 1.5, 'DisplayName', 'TTL OFF');
                else
                    xline(E.onsets(j), 'g-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
                    xline(E.falls(j),  'r-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
                end
                
                % Settle lines
                if j <= numel(E.settle_array) && ~isnan(E.settle_array(j))
                    settle_time_abs = E.falls(j) + E.settle_array(j);
                    if isempty(h_set) && i == 1
                        h_set = xline(settle_time_abs, 'b:', 'LineWidth', 2.0, ...
                                      'DisplayName', 'Settle Time');
                    else
                        xline(settle_time_abs, 'b:', 'LineWidth', 2.0, ...
                              'HandleVisibility', 'off');
                    end
                end
            end
            
            % Title with metrics
            title(sprintf('IBI: %.1f ms | %s | Rise: %.1f ms | Resid: %.1f%% | Final Settle: %.0f ms', ...
                E.ibi, E.Name, E.mean_rise, E.mean_resid, E.final_set), ...
                'Interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 10);
            ylabel('Acceleration (g)');
            if i == n_plots, xlabel('Time (ms) [0 = First Burst Onset]'); end
            
            % Dynamic limits
            x_end = E.falls(end) + max(300, E.final_set + 50);
            if isnan(x_end), x_end = E.falls(end) + 500; end
            xlim([-50, x_end]);
            yMax = max([max(abs(E.raw(:))), max(E.env(:))]);
            if yMax == 0 || isnan(yMax), yMax = 1; end
            ylim([-yMax*1.1, yMax*1.1]);
            
            % Legend on first subplot only
            if i == 1
                legs = [h_raw, h_env];
                if ~isempty(h_on),  legs = [legs, h_on];  end
                if ~isempty(h_off), legs = [legs, h_off]; end
                if ~isempty(h_set), legs = [legs, h_set]; end
                legend(legs, 'Location', 'northeast', 'Orientation', 'horizontal');
            end
        end
        sgtitle(sprintf('Train Dynamics: Pot %d — 250 Hz, 100ms bursts, 6 bursts/train', ...
                this_pot), 'FontSize', 13, 'FontWeight', 'bold');
    end
end