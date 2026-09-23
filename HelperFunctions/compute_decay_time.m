function decay_results = compute_decay_time(VibData, do_plot)
% COMPUTE_DECAY_TIME  Post-offset decay time constant and settling time.
%
%   decay_results = compute_decay_time(VibData)
%   decay_results = compute_decay_time(VibData, do_plot)
%
%   For each calibration burst, locates the stimulus-offset drop from the
%   known burst duration, fits a single exponential to the Hilbert envelope
%   over the 85%->35% decay band (log-linear fit), and measures settling
%   time to 5% of peak with an adaptive noise floor. Bursts below an SNR
%   gate are stored as NaN. Trains are skipped.
%
%   Output struct array fields:
%     .Position .PosIdx .PotStep .FreqHz — condition identifiers
%     .Tau_ms          — exponential decay time constant (ms)
%     .Settle_Meas_ms  — measured settling time to threshold (ms)
%     .R_Squared       — goodness of fit (log domain)
%     .BurstIdx .Source

    %% === PARAMETERS ===
    if nargin < 2, do_plot = true; end
    
    PLOT_PRE_MS      = 100;    
    PLOT_POST_MS     = 2000;    
    CALC_PADDING_MS  = 500;    
    
    FIT_UPPER_BOUND  = 0.85; 
    FIT_LOWER_BOUND  = 0.35; 
    BASE_SETTLE_TARGET = 0.05; 
    
    SNR_THRESHOLD_SIGMA = 6;
    MIN_FIT_MS          = 10;
    
    %% === PROCESS ALL BURSTS ===
    decay_results = struct('Position',{}, 'PosIdx',{}, 'PotStep',{}, ...
                           'FreqHz',{}, 'Tau_ms',{}, 'Settle_Meas_ms',{}, ...
                           'R_Squared',{}, 'BurstIdx',{}, 'Source',{});
    idx = 0;

    examples = struct();
    ex_count = 0;
    seen_combos = {};
    
    nPos = numel(VibData);
    fprintf('\n======================================================\n');
    fprintf('  DECAY ANALYSIS: KNOWN-DURATION DROP + SNR GATE\n');
    fprintf('======================================================\n');
    
    for p_idx = 1:nPos
        D = VibData(p_idx).Raw;
        if isempty(D), continue; end
        
        Fs = D.Fs;
        dt = 1/Fs;
        burst_list = get_burst_list(VibData(p_idx));
        
        % Global noise floor from pre-experiment baseline
        if numel(D.Z) > Fs
            global_noise_std = std(double(D.Z(1:round(0.1*Fs))));
        else
            global_noise_std = 0.001;
        end
        
        for bi = 1:numel(burst_list)
            B = burst_list{bi};
            if strcmp(B.source, 'Train'), continue; end
            
            % --- A. DROP INDEX FROM KNOWN BURST DURATION ---
            abs_drop_idx = B.onset + round(B.burst_ms / 1000 * Fs);
            if abs_drop_idx > numel(D.Z), continue; end
            
            % --- B. EXTRACT WINDOW ---
            preS  = round(PLOT_PRE_MS/1000 * Fs);
            postS = round(PLOT_POST_MS/1000 * Fs);
            padS  = round(CALC_PADDING_MS/1000 * Fs);
            
            w_start = abs_drop_idx - (preS + padS);
            w_end   = abs_drop_idx + (postS + padS);
            if w_start < 1 || w_end > numel(D.Z), continue; end
            
            wide_raw = double(D.Z(w_start:w_end));
            wide_raw = wide_raw - mean(wide_raw);
            try wide_env_raw = abs(hilbert(wide_raw)); catch, wide_env_raw = abs(wide_raw); end
            
            % --- C. FREQUENCY-DEPENDENT SMOOTHING ---
            if B.freq > 0
                smooth_win_ms = (5 / B.freq) * 1000;
            else
                smooth_win_ms = 20; 
            end
            smooth_win_ms = max(20, min(100, smooth_win_ms));
            wide_env = movmean(wide_env_raw, round(smooth_win_ms/1000 * Fs));
            
            % Crop to final view
            cut_i1 = padS + 1;
            cut_i2 = cut_i1 + preS + postS;
            final_raw = wide_raw(cut_i1:cut_i2);
            final_env = wide_env(cut_i1:cut_i2);
            
            t_axis = (-preS:postS) * dt * 1000; 
            t0_idx = preS + 1;
            
            % --- D. SNR GATE ---
            local_peak = final_env(t0_idx);
            if t0_idx + 50 <= numel(final_env) && t0_idx - 50 > 0
                peak_win = final_env(t0_idx-50:t0_idx+50);
                if local_peak < max(peak_win), local_peak = max(peak_win); end
            end
            
            if local_peak < SNR_THRESHOLD_SIGMA * global_noise_std
                idx = idx + 1;
                decay_results(idx).Position       = VibData(p_idx).Name;
                decay_results(idx).PosIdx         = p_idx;
                decay_results(idx).PotStep        = B.pot;
                decay_results(idx).FreqHz         = B.freq;
                decay_results(idx).Tau_ms         = NaN;
                decay_results(idx).Settle_Meas_ms = NaN;
                decay_results(idx).R_Squared      = NaN;
                decay_results(idx).BurstIdx       = B.bidx;
                decay_results(idx).Source         = B.source;
                continue;
            end
            
            % --- E. EXPONENTIAL FIT ---
            fit_upper = local_peak * FIT_UPPER_BOUND;
            fit_lower = max(local_peak * FIT_LOWER_BOUND, 3*global_noise_std);
            
            f_start = find(final_env(t0_idx:end) < fit_upper, 1, 'first');
            if isempty(f_start), f_start=0; end
            f_start = t0_idx + f_start - 1;
            
            f_end = find(final_env(f_start:end) < fit_lower, 1, 'first');
            if isempty(f_end)
                f_end = numel(final_env);
            else
                f_end = f_start + f_end - 1;
            end
            
            y_fit = final_env(f_start:f_end);
            t_fit = t_axis(f_start:f_end) / 1000;
            
            fit_span_ms = (f_end - f_start) * dt * 1000;
            min_fit_span = max(MIN_FIT_MS, 2000 / B.freq);
            
            tau_ms = NaN; r_sq = 0;
            if numel(y_fit) > 5 && local_peak > 1e-6 && fit_span_ms >= min_fit_span
                valid = y_fit > 1e-9;
                if sum(valid) > 5
                    p_poly = polyfit(t_fit(valid), log(y_fit(valid)), 1);
                    calculated_tau = -1000 / p_poly(1);
                    if calculated_tau > 2 && calculated_tau < 2000
                        tau_ms = calculated_tau;
                        y_pred = exp(p_poly(2) + p_poly(1)*t_fit(valid));
                        SSres = sum((log(y_fit(valid)) - log(y_pred)).^2);
                        SStot = sum((log(y_fit(valid)) - mean(log(y_fit(valid)))).^2);
                        r_sq = 1 - (SSres/SStot);
                    end
                end
            end
            
            % --- F. SETTLE TIME (protected tail window) ---
            tail_start_ms = 500;
            tail_end_ms   = min(800, PLOT_POST_MS);
            tail_i1 = t0_idx + round(tail_start_ms/1000*Fs);
            tail_i2 = t0_idx + round(tail_end_ms/1000*Fs);
            tail_i1 = min(tail_i1, numel(final_env));
            tail_i2 = min(tail_i2, numel(final_env));
            
            if tail_i2 > tail_i1
                tail_mean = mean(final_env(tail_i1:tail_i2));
                tail_std  = std(final_env(tail_i1:tail_i2));
            else
                tail_mean = 0; tail_std = global_noise_std;
            end
            
            adaptive_floor = tail_mean + 4 * tail_std;
            base_thresh    = local_peak * BASE_SETTLE_TARGET;
            settle_thresh  = max(base_thresh, adaptive_floor);
            
            settle_meas = NaN;
            dwell_samps = round(30/1000 * Fs);
            
            search_start = t0_idx + round(2/1000*Fs); 
            for k = search_start : (numel(final_env) - dwell_samps)
                if final_env(k) <= settle_thresh
                    if max(final_env(k:k+dwell_samps)) <= settle_thresh
                        settle_meas = t_axis(k);
                        break;
                    end
                end
            end
            
            % --- G. SETTLE SANITY CHECK ---
            if ~isnan(tau_ms) && ~isnan(settle_meas)
                if settle_meas < tau_ms * 0.5
                    settle_meas = NaN;
                end
            end
            
            % --- H. STORE ---
            idx = idx + 1;
            decay_results(idx).Position       = VibData(p_idx).Name;
            decay_results(idx).PosIdx         = p_idx;
            decay_results(idx).PotStep        = B.pot;
            decay_results(idx).FreqHz         = B.freq;
            decay_results(idx).Tau_ms         = tau_ms;
            decay_results(idx).Settle_Meas_ms = settle_meas;
            decay_results(idx).R_Squared      = r_sq;
            decay_results(idx).BurstIdx       = B.bidx;
            decay_results(idx).Source         = B.source;
            
            % Save example
            combo_key = sprintf('P%d_F%d', B.pot, B.freq);
            if ~ismember(combo_key, seen_combos)
                seen_combos{end+1} = combo_key;
                ex_count = ex_count + 1;
                examples(ex_count).t          = t_axis;
                examples(ex_count).raw        = final_raw;
                examples(ex_count).env        = final_env; 
                examples(ex_count).thresh     = settle_thresh;
                examples(ex_count).tau        = tau_ms;
                examples(ex_count).settle_meas = settle_meas;
                examples(ex_count).pot        = B.pot;
                examples(ex_count).freq       = B.freq;
                examples(ex_count).source     = B.source;
                examples(ex_count).peak_amp   = local_peak;
            end
        end
        fprintf('  %s: done\n', VibData(p_idx).Name);
    end
    
    %% === REPORTING (always printed) ===
    if idx == 0, warning('No decay times computed.'); return; end
    
    fprintf('\n=================================================================================\n');
    fprintf('  DETAILED BURST ANALYSIS LOG\n');
    fprintf('=================================================================================\n');
    fprintf('%-12s | %-6s | %-4s | %-4s | %-10s | %-10s\n', ...
            'SOURCE', 'FREQ', 'POT', 'IDX', 'TAU (ms)', 'SETTLE (ms)');
    fprintf('%s\n', repmat('-', 1, 65));
    for i = 1:idx
        res = decay_results(i);
        fprintf('%-12s | %3d Hz | P%-3d | #%-3d | %8.2f   | %8.2f\n', ...
            res.Source, res.FreqHz, res.PotStep, res.BurstIdx, res.Tau_ms, res.Settle_Meas_ms);
    end
    fprintf('%s\n', repmat('-', 1, 65));

    all_tau  = [decay_results.Tau_ms];
    all_set  = [decay_results.Settle_Meas_ms];
    all_freq = [decay_results.FreqHz];
    all_pot  = [decay_results.PotStep];
    
    fprintf('\n=========================================================================\n');
    fprintf('  TABLE 1: DECAY SUMMARY BY FREQUENCY\n');
    fprintf('=========================================================================\n');
    fprintf('%-8s | %-15s | %-15s | %-6s\n', 'FREQ', 'Tau (ms)', 'Settle (ms)', 'N');
    fprintf('%s\n', repmat('-', 1, 52));
    u_freq = unique(all_freq);
    for j = 1:numel(u_freq)
        f = u_freq(j); mask = (all_freq == f);
        dTau = all_tau(mask); dSet = all_set(mask);
        fprintf('%3d Hz   | %6.2f +/- %-6.1f | %6.2f +/- %-6.1f | %d\n', ...
            f, mean(dTau,'omitnan'), std(dTau,'omitnan'), ...
               mean(dSet,'omitnan'), std(dSet,'omitnan'), sum(mask));
    end
    
    fprintf('\n=========================================================================\n');
    fprintf('  TABLE 2: DECAY BY FREQUENCY x POT STEP\n');
    fprintf('=========================================================================\n');
    fprintf('%-8s | %-6s | %-15s | %-15s | %-6s\n', 'FREQ', 'POT', 'Tau (ms)', 'Settle (ms)', 'N');
    fprintf('%s\n', repmat('-', 1, 62));
    u_pot = unique(all_pot);
    for j = 1:numel(u_freq)
        for pi = 1:numel(u_pot)
            f = u_freq(j); p = u_pot(pi);
            mask = (all_freq == f) & (all_pot == p);
            if sum(mask) == 0, continue; end
            dTau = all_tau(mask); dSet = all_set(mask);
            fprintf('%3d Hz   | P%-4d | %6.2f +/- %-6.1f | %6.2f +/- %-6.1f | %d\n', ...
                f, p, mean(dTau,'omitnan'), std(dTau,'omitnan'), ...
                      mean(dSet,'omitnan'), std(dSet,'omitnan'), sum(mask));
        end
    end
    
    %% === PLOTS (conditional) ===
    if do_plot && ex_count > 0
        plot_decay_diagnostic(examples, decay_results);
    end
end

%% ========================================================================
function plot_decay_diagnostic(ex, all_res)
    full_freqs = [all_res.FreqHz];
    full_pots  = [all_res.PotStep];
    full_set   = [all_res.Settle_Meas_ms];
    
    ex_freqs = [ex.freq];
    unique_freqs = sort(unique(ex_freqs), 'descend');
    
    for fi = 1:numel(unique_freqs)
        this_freq = unique_freqs(fi);
        mask = (ex_freqs == this_freq);
        curr_ex = ex(mask);
        if isempty(curr_ex), continue; end
        
        [~, sort_idx] = sort([curr_ex.pot]);
        curr_ex = curr_ex(sort_idx);
        
        num_plots = min(7, numel(curr_ex));
        ncols = min(4, num_plots);
        nrows = ceil(num_plots / ncols);
        
        figure('Name', sprintf('Decay: %d Hz', this_freq), ...
               'Color', 'w', 'Position', [60+fi*40, 60+fi*40, 350*ncols, 350*nrows]);
        
        for i = 1:num_plots
            subplot(nrows, ncols, i); hold on; grid on; box on;
            E = curr_ex(i);
            
            cond_mask = (full_freqs == E.freq) & (full_pots == E.pot);
            mean_settle = mean(full_set(cond_mask), 'omitnan');
            
            plot(E.t, E.raw, 'Color', [0.85 0.85 0.85], 'LineWidth', 0.5);
            plot(E.t, E.env, 'k-', 'LineWidth', 1.2); 
            xline(0, 'k-', 'LineWidth', 2.0);
            if ~isnan(E.settle_meas), xline(E.settle_meas, 'r-', 'LineWidth', 1.5); end
            if ~isnan(mean_settle), xline(mean_settle, 'r:', 'LineWidth', 2.0); end
            yline(E.thresh, 'b:', 'LineWidth', 0.5);
            
            info_str = {
                sprintf('\\bfPot %d\\rm', E.pot);
                sprintf('Tau: %.1f ms', E.tau);
                sprintf('Settle: \\color{red}%.0f ms', E.settle_meas);
                sprintf('Mean: \\color{red}%.0f ms', mean_settle);
            };
            text(0.55, 0.78, info_str, 'Units', 'normalized', ...
                'BackgroundColor', [1 1 1 0.8], 'EdgeColor', 'k', 'FontSize', 9);
            
            xlim([E.t(1), E.t(end)]);
            yMax = max([max(abs(E.raw(:))), max(E.env(:))]);
            if yMax > 0, ylim([-yMax*1.1, yMax*1.1]); end
            
            xlabel('Time (ms) [0 = Burst End]');
            if mod(i-1, ncols) == 0, ylabel('Acceleration (g)'); end
            title(sprintf('Pot %d', E.pot), 'FontSize', 11);
        end
        sgtitle(sprintf('Decay Analysis: %d Hz  (Tau fit 85%%->35%%, Settle to 5%%)', ...
                this_freq), 'FontSize', 13, 'FontWeight', 'bold');
    end
end

%% ========================================================================
function list = get_burst_list(S)
    list = {};
    if isfield(S, 'Calibration') && ~isempty(S.Calibration)
        for c = 1:numel(S.Calibration)
            C = S.Calibration(c);
            for k = 1:C.NBursts
                list{end+1} = struct('onset', C.OnsetIdx(k), 'pot', C.PotStep, ...
                    'freq', C.FreqHz, 'bidx', k, 'source', 'Calibration', ...
                    'burst_ms', C.BurstMs); %#ok<AGROW>
            end
        end
    end
    if isfield(S, 'Trains') && ~isempty(S.Trains)
        for c = 1:numel(S.Trains)
            T = S.Trains(c);
            if T.NBursts >= 1
                list{end+1} = struct('onset', T.BurstOnsetIdx(1), 'pot', T.PotStep, ...
                    'freq', T.FreqHz, 'bidx', 1, 'source', 'Train', ...
                    'burst_ms', T.BurstMs); %#ok<AGROW>
            end
        end
    end
end