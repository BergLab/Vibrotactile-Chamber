function rise_results = compute_rise_time(VibData, do_plot)
% COMPUTE_RISE_TIME  Envelope rise time of vibrotactile calibration bursts.
%
%   rise_results = compute_rise_time(VibData)
%   rise_results = compute_rise_time(VibData, do_plot)
%
%   For each calibration burst, detects the mechanical onset in the Z-axis
%   acceleration, extracts the Hilbert envelope of the low-pass filtered
%   signal, and measures rise time to 90% of the steady-state amplitude
%   (median envelope over the late 60–95% of the burst). Two rise times are
%   reported: from mechanical onset, and from the 10% crossing.
%
%   Output struct array fields:
%     .Position .PotStep .FreqHz   — condition identifiers
%     .RiseTime_Mech90_ms  — mechanical onset to 90% (ms)
%     .RiseTime_10_90_ms   — 10% to 90% rise time (ms)
%     .BurstIdx
    %% === PARAMETERS ===
    if nargin < 2, do_plot = true; end
    PRE_MS        = 20;
    POST_MS       = 300;
    
    LOW_SIGMA     = 2.0;    
    HIGH_SIGMA    = 8.0; 
    SMOOTH_CYCLES = 1.0;    

    %% === INITIALIZATION ===
    rise_results = struct('Position',{}, 'PotStep',{}, 'FreqHz',{}, ...
                          'RiseTime_Mech90_ms',{}, 'RiseTime_10_90_ms',{}, ...
                          'BurstIdx',{});
    idx = 0;
    
    examples = struct();
    ex_count = 0;
    seen_combos = {}; 
    
    nPos = numel(VibData);
    fprintf('Computing Rise Times (Envelope Analysis)...\n');

    for p = 1:nPos
        D = VibData(p).Raw;
        if isempty(D), continue; end
        
        Fs = D.Fs;
        dt = 1/Fs;
        pre_samp = round(PRE_MS/1000 * Fs);
        
        burst_list = get_burst_list(VibData(p)); 
        
        for bi = 1:numel(burst_list)
            B = burst_list{bi};
            
            % --- 1. Extract Window ---
            i1 = B.onset - pre_samp;
            i2 = B.onset + round(POST_MS/1000 * Fs);
            if i1 < 1 || i2 > numel(D.Z), continue; end
            
            raw_seg = double(D.Z(i1:i2));
            t_axis = ((0:numel(raw_seg)-1) - pre_samp) * dt * 1000;
            
            % --- 2. Pre-processing ---
            bl_indices = 1:pre_samp;
            bl_mean = mean(raw_seg(bl_indices));
            bl_std  = std(raw_seg(bl_indices));
            centered_seg = raw_seg - bl_mean;
            
            % --- 3. Envelope Extraction ---
            [b_lp, a_lp] = butter(2, 500/(Fs/2), 'low');
            clean_sig = filtfilt(b_lp, a_lp, centered_seg);
            
            analytic_sig = hilbert(clean_sig);
            env_raw = abs(analytic_sig);
            period_samps = round(Fs / B.freq * SMOOTH_CYCLES);
            envelope = movmean(env_raw, period_samps);
            
            % --- 4. Find Mechanical Onset (t_mech) ---
            low_thresh = LOW_SIGMA * bl_std;
            high_thresh = HIGH_SIGMA * bl_std;
            
            mech_idx = [];
            search_end = pre_samp + round(30/1000*Fs); 
            
            for k = (pre_samp+1):search_end
                if abs(centered_seg(k)) > low_thresh
                    future_peak = max(abs(centered_seg(k:min(end, k+round(10/1000*Fs)))));
                    if future_peak > high_thresh
                        for b = k:-1:1
                            if abs(centered_seg(b)) < 1.5*bl_std
                                mech_idx = b;
                                break;
                            end
                        end
                        if isempty(mech_idx), mech_idx = k; end
                        break; 
                    end
                end
            end
            
            if isempty(mech_idx), continue; end 
            
            % --- 5. Determine Steady State Amplitude (100%) ---
            % Use identical 60% to 95% late window logic
            burst_samp = round(B.burst_ms / 1000 * Fs);
            ss_start = pre_samp + round(burst_samp * 0.60);
            ss_end   = pre_samp + round(burst_samp * 0.95);
            
            % Failsafe if the extraction window (POST_MS) is shorter than the burst
            if ss_end > numel(envelope)
                ss_end = numel(envelope);
            end
            
            if ss_start >= ss_end
                target_amp = max(envelope);
            else
                target_amp = median(envelope(ss_start:ss_end));
            end
            
            % --- 6. Calculate Rise Times ---
            threshold_10 = 0.10 * target_amp;
            threshold_90 = 0.90 * target_amp;
            
            idx_10 = find(envelope(mech_idx:end) >= threshold_10, 1, 'first');
            if isempty(idx_10), continue; end
            idx_10 = idx_10 + mech_idx - 1;
            
            idx_90 = find(envelope(mech_idx:end) >= threshold_90, 1, 'first');
            if isempty(idx_90), continue; end
            idx_90 = idx_90 + mech_idx - 1;
            
            t_mech = t_axis(mech_idx);
            t_90   = t_axis(idx_90);
            t_10   = t_axis(idx_10);
            
            rt_mech_90 = t_90 - t_mech;
            rt_10_90 = t_90 - t_10;
            
            % --- 7. Store Results ---
            idx = idx + 1;
            rise_results(idx).Position = VibData(p).Name;
            rise_results(idx).PotStep  = B.pot;
            rise_results(idx).FreqHz   = B.freq;
            rise_results(idx).RiseTime_Mech90_ms = rt_mech_90;
            rise_results(idx).RiseTime_10_90_ms  = rt_10_90;
            rise_results(idx).BurstIdx = B.bidx;
            
            % --- 8. Store Example ---
            combo_key = sprintf('%d_%d', B.pot, B.freq);
            if ~ismember(combo_key, seen_combos)
                seen_combos{end+1} = combo_key;
                ex_count = ex_count + 1;
                examples(ex_count).t = t_axis;
                examples(ex_count).raw = centered_seg;
                examples(ex_count).t_mech = t_mech;
                examples(ex_count).t_90 = t_90;
                examples(ex_count).pot = B.pot;
                examples(ex_count).freq = B.freq;
            end
        end
    end
    
    fprintf('Done. Processed %d bursts.\n', idx);
    
    if idx > 0
        report_rise_stats(rise_results);
    end
    
    if do_plot && ex_count > 0
        plot_rise_diagnostic(examples);
    end
end

function report_rise_stats(res)
    all_pot = [res.PotStep];
    all_freq = [res.FreqHz];
    all_mech90 = [res.RiseTime_Mech90_ms];
    
    fprintf('\n===== RISE TIME (Mech Onset -> 90%%) =====\n');
    fprintf('%-10s %-10s %-12s %-12s %-10s\n', 'Pot Step', 'Freq', 'Mean (ms)', 'SD (ms)', 'Count');
    fprintf('%s\n', repmat('-', 1, 56));
    
    combos = unique([all_pot', all_freq'], 'rows');
    for i = 1:size(combos, 1)
        pot = combos(i,1); freq = combos(i,2);
        mask = (all_pot == pot) & (all_freq == freq);
        vals = all_mech90(mask);
        fprintf('%-10d %-10d %-12.3f %-12.3f %-10d\n', pot, freq, mean(vals), std(vals), numel(vals));
    end
end

function plot_rise_diagnostic(ex)
    figure('Name', 'Rise Time Check', 'Color', 'w', 'Position', [100 100 1400 800]);
    num_plots = min(6, numel(ex));
    
    for i = 1:num_plots
        subplot(2,3,i); hold on; box on;
        E = ex(i);
        
        plot(E.t, E.raw, 'k-', 'LineWidth', 0.8);
        xline(0, '--k', 'LineWidth', 1.0);
        xline(E.t_mech, '-g', 'LineWidth', 2.0);
        xline(E.t_90, '-r', 'LineWidth', 2.0);
        
        rise_val = E.t_90 - E.t_mech;
        title(sprintf('Pot %d | %d Hz\nRise: %.1f ms', E.pot, E.freq, rise_val));
        xlabel('Time (ms) rel. to TTL');
        ylabel('Acceleration (g)');
        
        xlim([-10, 60]); 
        mask = E.t >= -10 & E.t <= 60;
        if any(mask)
            ymax = max(abs(E.raw(mask)));
            ylim([-ymax*1.1, ymax*1.1]);
        end
    end
end

% === HELPER: GET BURST LIST (v3: Calibration only) ===
function list = get_burst_list(S)
    list = {};
    if isfield(S, 'Calibration') && ~isempty(S.Calibration)
        for c = 1:numel(S.Calibration)
            C = S.Calibration(c);
            for k = 1:C.NBursts
                % ADDED 'burst_ms' to the struct here:
                list{end+1} = struct('onset', C.OnsetIdx(k), 'pot', C.PotStep, ...
                    'freq', C.FreqHz, 'bidx', k, 'burst_ms', C.BurstMs); 
            end
        end
    end
end