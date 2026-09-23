function delay_results = compute_mechanical_delay(VibData, do_plot)
% COMPUTE_MECHANICAL_DELAY  Latency from stimulus onset (TTL) to mechanical onset.
%
%   delay_results = compute_mechanical_delay(VibData)
%   delay_results = compute_mechanical_delay(VibData, do_plot)
%
%   For each burst (calibration and the first burst of each train), detects
%   the mechanical onset in the Z-axis acceleration relative to the commanded
%   TTL onset. A two-stage detector finds a trigger (amplitude crossing
%   validated by a downstream peak), then backtracks along amplitude and slope
%   to the true onset; a minimum physical delay gates false-early detections.
%
%   Output struct array fields:
%     .Position .PosIdx .PotStep .FreqHz — condition identifiers
%     .Delay_ms  — onset latency relative to TTL (ms)
%     .BurstIdx .Source
    %% === PARAMETERS ===
    if nargin < 2, do_plot = true; end
    PRE_MS       = 50;    
    POST_MS      = 100;   
    
    % --- PHYSICAL CONSTRAINTS ---
    MIN_MECH_DELAY_MS = 0.5; 
    
    % --- DETECTION SETTINGS ---
    LOW_SIGMA    = 2.0;   
    HIGH_SIGMA   = 6.0;   
    VALIDATION_MS = 15;   
    
    % --- BACKTRACK SETTINGS ---
    AMP_FLOOR_SIGMA   = 1.5; 
    SLOPE_FLOOR_SIGMA = 2.0; 
    
    %% === PROCESS ALL BURSTS ===
    delay_results = struct('Position',{}, 'PosIdx',{}, 'PotStep',{}, ...
                           'FreqHz',{}, 'Delay_ms',{}, 'BurstIdx',{}, ...
                           'Source',{});
    idx = 0;
    
    examples = struct();
    ex_count = 0;
    seen_combos = {};
    
    nPos = numel(VibData);
    fprintf('\n======================================================\n');
    fprintf('  MECHANICAL DELAY ANALYSIS (Robust Low-Freq Fix)\n');
    fprintf('======================================================\n');
    fprintf('Processing bursts across %d positions...\n', nPos);
    
    for p = 1:nPos
        D = VibData(p).Raw;
        if isempty(D), continue; end
        
        Fs = D.Fs;
        dt = 1/Fs;
        pre_samp = round(PRE_MS/1000 * Fs);
        lookahead_samps = round(VALIDATION_MS/1000 * Fs);
        min_delay_samps = round(MIN_MECH_DELAY_MS/1000 * Fs);
        
        burst_list = get_burst_list(VibData(p)); 
        
        for bi = 1:numel(burst_list)
            B = burst_list{bi};
            
            % --- 1. Extract Window ---
            i1 = B.onset - pre_samp;
            i2 = B.onset + round(POST_MS/1000 * Fs);
            if i1 < 1 || i2 > numel(D.Z), continue; end
            
            raw_seg = double(D.Z(i1:i2));
            t_axis = ((0:numel(raw_seg)-1) - pre_samp) * dt * 1000; 
            
            % --- 2. Baseline Statistics ---
            bl_indices = 1:(pre_samp - round(2/1000*Fs)); 
            if isempty(bl_indices), bl_indices = 1:pre_samp; end
            
            bl_mean = mean(raw_seg(bl_indices));
            bl_std  = std(raw_seg(bl_indices));
            if bl_std < 1e-9, bl_std = 1e-9; end 
            
            centered_seg = raw_seg - bl_mean;
            
            slope_sig = [0; diff(centered_seg)]; 
            slope_std = std(slope_sig(bl_indices));
            
            low_thresh  = LOW_SIGMA * bl_std;
            high_thresh = HIGH_SIGMA * bl_std;
            amp_floor   = AMP_FLOOR_SIGMA * bl_std;
            
            if B.freq < 100
                slope_scale = 0.5; 
            else
                slope_scale = 1.0;
            end
            slope_floor = (SLOPE_FLOOR_SIGMA * slope_scale) * slope_std;
            
            % --- 3. Find Trigger (Scout) ---
            search_start = pre_samp + min_delay_samps + 1; 
            search_end   = min(numel(centered_seg) - lookahead_samps, ...
                               pre_samp + round(40/1000*Fs)); 
            
            trigger_idx = [];
            for k = search_start:search_end
                if abs(centered_seg(k)) > low_thresh
                    future_peak = max(abs(centered_seg(k:min(end, k+lookahead_samps))));
                    if future_peak > high_thresh
                        trigger_idx = k;
                        break; 
                    end
                end
            end
            
            if isempty(trigger_idx), continue; end
            
            % --- 4. Precision Backtrack (Slope) ---
            onset_idx = trigger_idx;
            max_backtrack_ms = (1000 / B.freq) * 0.6; 
            
            for k = trigger_idx:-1:1
                if t_axis(k) <= MIN_MECH_DELAY_MS
                    onset_idx = k; 
                    break;
                end

                val_amp = abs(centered_seg(k));
                val_slope = abs(slope_sig(k));
                
                if val_amp <= amp_floor || val_slope <= slope_floor
                    onset_idx = k; 
                    break;
                end
                
                if (t_axis(trigger_idx) - t_axis(k)) > max_backtrack_ms
                    onset_idx = k; 
                    break;
                end
            end
            
            delay_ms = t_axis(onset_idx);
            
            % --- 5. SANITY CHECK CORRECTION ---
            if delay_ms <= MIN_MECH_DELAY_MS + 0.05
                start_search_idx = pre_samp + min_delay_samps;
                robust_thresh = 3.0 * bl_std;
                rel_idx = find(abs(centered_seg(start_search_idx:trigger_idx)) > robust_thresh, 1);
                
                if ~isempty(rel_idx)
                    onset_idx = start_search_idx + rel_idx - 1;
                    delay_ms  = t_axis(onset_idx);
                else
                    delay_ms = MIN_MECH_DELAY_MS;
                end
            end
            
            % --- 6. Store Result ---
            idx = idx + 1;
            delay_results(idx).Position = VibData(p).Name;
            delay_results(idx).PosIdx   = p;
            delay_results(idx).PotStep  = B.pot;
            delay_results(idx).FreqHz   = B.freq;
            delay_results(idx).Delay_ms = delay_ms;
            delay_results(idx).BurstIdx = B.bidx;
            delay_results(idx).Source   = B.source;
            
            % --- 7. Store Plotting Example ---
            combo_key = sprintf('%s_P%d_F%d', B.source, B.pot, B.freq);
            if ~ismember(combo_key, seen_combos)
                seen_combos{end+1} = combo_key; 
                ex_count = ex_count + 1;
                examples(ex_count).t = t_axis;
                examples(ex_count).raw = raw_seg;
                examples(ex_count).bl_mean = bl_mean;
                examples(ex_count).amp_floor = amp_floor;
                examples(ex_count).onset_t = delay_ms;
                examples(ex_count).onset_idx = onset_idx;
                examples(ex_count).pot = B.pot;
                examples(ex_count).freq = B.freq;
                examples(ex_count).source = B.source;
            end
        end
        fprintf('  %s: done\n', VibData(p).Name);
    end
    
    fprintf('\nProcessed %d total bursts.\n', idx);
    %% === STATISTICS & REPORTING ===
    if idx == 0, warning('No delays computed.'); return; end
    
    all_pos   = {delay_results.Position};
    all_pidx  = [delay_results.PosIdx];
    all_pot   = [delay_results.PotStep];
    all_freq  = [delay_results.FreqHz];
    all_delay = [delay_results.Delay_ms];
    all_src   = {delay_results.Source};
    
    print_tables(all_src, all_delay, all_pot, all_freq, all_pidx, all_pos, delay_results);
    
    if do_plot && ex_count > 0
        plot_diagnostic_extended(examples);
    end
end

%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================
function print_tables(all_src, all_delay, all_pot, all_freq, all_pidx, ~, delay_results)
    % TABLE 1: GRAND SUMMARY
    fprintf('\n======================================================\n');
    fprintf('  TABLE 1: GRAND SUMMARY BY SOURCE\n');
    fprintf('======================================================\n');
    fprintf('%-12s %-10s %-10s %-10s %-8s %-8s\n', 'Source', 'Mean', 'SD', 'Med', 'IQR', 'N');
    fprintf('%s\n', repmat('-', 1, 60));
    u_src = unique(all_src);
    for i = 1:numel(u_src)
        s = u_src{i}; mask = strcmp(all_src, s); d = all_delay(mask);
        fprintf('%-12s %-10.3f %-10.3f %-10.3f %-8.3f %-8d\n', s, mean(d), std(d), median(d), iqr(d), numel(d));
    end
    
    % TABLE 2: CALIBRATION BY POT x FREQ
    mask_cal = strcmp(all_src, 'Calibration');
    if any(mask_cal)
        fprintf('\n======================================================\n');
        fprintf('  TABLE 2: CALIBRATION BY POT x FREQ\n');
        fprintf('======================================================\n');
        fprintf('%-10s %-10s %-10s %-10s %-10s %-8s\n', 'Pot', 'Freq', 'Mean', 'SD', 'Med', 'N');
        fprintf('%s\n', repmat('-', 1, 68));
        cal_combos = unique([all_pot(mask_cal)', all_freq(mask_cal)'], 'rows');
        for i = 1:size(cal_combos, 1)
            p = cal_combos(i,1); f = cal_combos(i,2);
            m = mask_cal & (all_pot == p) & (all_freq == f); d = all_delay(m);
            fprintf('%-10d %-10d %-10.3f %-10.3f %-10.3f %-8d\n', p, f, mean(d), std(d), median(d), numel(d));
        end
    end
    
    % TABLE 3: TRAINS
    mask_tr = strcmp(all_src, 'Train');
    if any(mask_tr)
        fprintf('\n======================================================\n');
        fprintf('  TABLE 3: TRAINS BY POT\n');
        fprintf('======================================================\n');
        fprintf('%-10s %-10s %-10s %-10s %-8s\n', 'Pot', 'Mean', 'SD', 'Med', 'N');
        fprintf('%s\n', repmat('-', 1, 48));
        u_pot = unique(all_pot(mask_tr));
        for i = 1:numel(u_pot)
            m = mask_tr & (all_pot == u_pot(i)); d = all_delay(m);
            fprintf('%-10d %-10.3f %-10.3f %-10.3f %-8d\n', u_pot(i), mean(d), std(d), median(d), numel(d));
        end
    end
    
    % TABLE 4: PER POSITION
    fprintf('\n======================================================\n');
    fprintf('  TABLE 4: PER-POSITION SUMMARY\n');
    fprintf('======================================================\n');
    fprintf('%-10s %-10s %-10s %-10s %-8s\n', 'Position', 'Mean', 'SD', 'Med', 'N');
    fprintf('%s\n', repmat('-', 1, 58));
    u_p = unique(all_pidx);
    for i = 1:numel(u_p)
        mask = (all_pidx == u_p(i)); d = all_delay(mask);
        pname = delay_results(find(mask,1)).Position;
        fprintf('%-10s %-10.3f %-10.3f %-10.3f %-8d\n', pname, mean(d), std(d), median(d), numel(d));
    end
end

function plot_diagnostic_extended(ex)
    all_sources = {ex.source};
    unique_sources = unique(all_sources);
    
    for s = 1:numel(unique_sources)
        src_name = unique_sources{s};
        mask = strcmp(all_sources, src_name);
        curr_ex = ex(mask);
        [~, sort_idx] = sortrows([[curr_ex.freq]', [curr_ex.pot]']);
        curr_ex = curr_ex(sort_idx);
        
        figure('Name', sprintf('Delay: %s', src_name), 'Color', 'w', 'Position', [50+(s*30) 50+(s*30) 1200 700]);
        num_plots = min(6, numel(curr_ex));
        
        for i = 1:num_plots
            subplot(2, 3, i); hold on; grid on;
            E = curr_ex(i);
            
            if E.freq <= 100
                x_lims = [-5, 20];
            elseif E.freq <= 150
                x_lims = [-2, 12];
            else
                x_lims = [-1, 5];
            end
            
            plot(E.t, E.raw, 'k-', 'LineWidth', 0.8);
            xline(0, 'k--', 'LineWidth', 1.0); 
            xline(E.onset_t, 'g-', 'LineWidth', 2.0);
            yline(E.bl_mean + E.amp_floor, ':', 'Color', [0.5 0.5 0.5]);
            
            xlim(x_lims);
            
            mask_t = (E.t >= x_lims(1)) & (E.t <= x_lims(2));
            local_data = E.raw(mask_t);
            if ~isempty(local_data)
                y_min = min(local_data); y_max = max(local_data);
                margin = max(0.1, (y_max - y_min) * 0.2);
                ylim([y_min - margin, y_max + margin]);
            end
            
            title(sprintf('%s | Pot %d | %dHz\nDelay: %.3f ms', E.source, E.pot, E.freq, E.onset_t), 'Interpreter', 'none');
        end
        sgtitle(sprintf('Source: %s', src_name));
    end
end

function list = get_burst_list(S)
    list = {};
    % Calibration bursts (replaces Amplitude + FreqSweep)
    if isfield(S, 'Calibration') && ~isempty(S.Calibration)
        for c = 1:numel(S.Calibration)
            C = S.Calibration(c);
            for k = 1:C.NBursts
                list{end+1} = struct('onset', C.OnsetIdx(k), 'pot', C.PotStep, ...
                    'freq', C.FreqHz, 'bidx', k, 'source', 'Calibration'); %#ok<AGROW>
            end
        end
    end
    % Train bursts (first burst of each train)
    if isfield(S, 'Trains') && ~isempty(S.Trains)
        for c = 1:numel(S.Trains)
            T = S.Trains(c);
            if T.NBursts >= 1
                list{end+1} = struct('onset', T.BurstOnsetIdx(1), 'pot', T.PotStep, ...
                    'freq', T.FreqHz, 'bidx', 1, 'source', 'Train'); %#ok<AGROW>
            end
        end
    end
end