function disp_results = compute_displacement(VibData, do_plot)
% COMPUTE_DISPLACEMENT  Compute calibrated displacement (um) per burst.
%
%   disp_results = compute_displacement(VibData)
%   disp_results = compute_displacement(VibData, do_plot)
%
%   For each calibration burst:
%     1. Double-integrates Z-axis acceleration to displacement
%     2. Computes Hilbert envelope on wide padded window (no edge artifacts)
%     3. Defines a steady-state window over the late 60–95% of the burst,
%        where mechanical transients have settled
%     4. Reports RMS and peak-to-peak displacement over the plateau
%
%   Output struct array fields:
%     .Position          — 'Pos1'..'Pos5'
%     .PosIdx            — position index
%     .PotStep           — pot step value
%     .FreqHz            — stimulus frequency
%     .BurstIdx          — rep number within condition
%     .Displacement_um   — peak-to-peak displacement (um)
%     .RMS_um            — RMS displacement (um)
%     .RMS_g             — RMS acceleration (g, middle 80%)
%     .SteadyStart_ms    — 10–90% rise time of the displacement envelope (ms)
%     .SteadyDur_ms      — duration of steady-state window (ms)
%     .Source            — 'Calibration'

    if nargin < 2, do_plot = true; end

    %% === PARAMETERS ===
    HPF_CUTOFF = 20;       % Hz — high-pass for integration drift removal
    HPF_ORDER  = 2;        % filter order (4th with filtfilt)
    G_TO_MS2   = 9.81;
    WIDE_PAD_MS = 500;     % padding to absorb filter/Hilbert edge effects

    %% === PROCESS ALL BURSTS ===
    disp_results = struct('Position',{}, 'PosIdx',{}, 'PotStep',{}, ...
                          'FreqHz',{}, 'BurstIdx',{}, ...
                          'Displacement_um',{}, 'RMS_um',{}, 'RMS_g',{}, ...
                          'SteadyStart_ms',{}, 'SteadyDur_ms',{}, 'Source',{});
    idx = 0;
    nPos = numel(VibData);

    fprintf('\n======================================================\n');
    fprintf('  DISPLACEMENT ANALYSIS: Rise-100%% to TTL-off\n');
    fprintf('======================================================\n');

    for p = 1:nPos
        D = VibData(p).Raw;
        if isempty(D), continue; end

        Fs = D.Fs;
        dt = 1 / Fs;
        [b_hp, a_hp] = butter(HPF_ORDER, HPF_CUTOFF / (Fs/2), 'high');

        if ~isfield(VibData(p), 'Calibration') || isempty(VibData(p).Calibration)
            continue;
        end

        for c = 1:numel(VibData(p).Calibration)
            C = VibData(p).Calibration(c);
            for k = 1:C.NBursts
                onset_idx  = C.OnsetIdx(k);
                burst_ms   = C.BurstMs;
                burst_samp = round(burst_ms / 1000 * Fs);
                offset_idx = onset_idx + burst_samp;

                % --- Wide window extraction (with padding) ---
                pad_samp = round(WIDE_PAD_MS / 1000 * Fs);
                i1 = max(1, onset_idx - pad_samp);
                i2 = min(numel(D.Z), offset_idx + pad_samp);
                if i2 - i1 < burst_samp * 0.5, continue; end

                seg_wide = double(D.Z(i1:i2));
                seg_wide = seg_wide - mean(seg_wide);

                % --- Acceleration RMS (middle 80% of raw burst) ---
                local_onset  = onset_idx - i1 + 1;
                local_offset = offset_idx - i1 + 1;
                margin = round(burst_samp * 0.1);
                acc_i1 = local_onset + margin;
                acc_i2 = local_offset - margin;

                if acc_i2 > acc_i1
                    acc_steady = seg_wide(acc_i1:acc_i2);
                    acc_steady = acc_steady - mean(acc_steady);
                    rms_g = sqrt(mean(acc_steady.^2));
                else
                    rms_g = NaN;
                end

                % --- Double integration on wide segment ---
                seg_filt = filtfilt(b_hp, a_hp, seg_wide);
                accel_ms2 = seg_filt * G_TO_MS2;

                velocity = cumtrapz(accel_ms2) * dt;
                velocity = filtfilt(b_hp, a_hp, velocity);

                displacement_m = cumtrapz(velocity) * dt;
                displacement_m = filtfilt(b_hp, a_hp, displacement_m);

                disp_um = displacement_m * 1e6;

                % Baseline-zero using 200ms pre-onset window
                bl_start = max(1, local_onset - round(200/1000*Fs));
                bl_end   = max(1, local_onset - round(10/1000*Fs));
                if bl_end > bl_start
                    disp_um = disp_um - mean(disp_um(bl_start:bl_end));
                end

                % --- Envelope on wide window (Hilbert edge-free) ---
                env_raw = abs(hilbert(disp_um));
                smooth_samps = max(round(20/1000*Fs), round(5/C.FreqHz*Fs));
                env = movmean(env_raw, smooth_samps);

                % --- Find steady-state amplitude (plateau) ---
                % Use the late window (60% to 95% of burst) to ensure mechanical transients have completely settled
                late_start = local_onset + round(burst_samp * 0.60);
                late_end   = local_onset + round(burst_samp * 0.95);
                
                if late_end <= late_start
                    late_start = local_onset + margin;
                    late_end   = local_offset - margin;
                end
                
                if late_end > late_start
                    plateau_amp = median(env(late_start:late_end));
                else
                    plateau_amp = max(env(local_onset:local_offset));
                end

                % --- Calculate True 10-90% Rise Time ---
                thresh_10 = 0.10 * plateau_amp;
                thresh_90 = 0.90 * plateau_amp;
                
                idx_10_local = find(env(local_onset:local_offset) >= thresh_10, 1, 'first');
                idx_90_local = find(env(local_onset:local_offset) >= thresh_90, 1, 'first');
                
                if ~isempty(idx_10_local) && ~isempty(idx_90_local) && (idx_90_local > idx_10_local)
                    steady_start_ms = (idx_90_local - idx_10_local) * dt * 1000; 
                else
                    steady_start_ms = NaN;
                end
                
                % Duration of the guaranteed steady-state window
                steady_dur_ms = (late_end - late_start) * dt * 1000;

                % --- Displacement over GUARANTEED steady-state window ---
                % We safely use the late window (60% to 95%) defined earlier for the plateau
                ss_i1 = late_start;
                ss_i2 = late_end;
                
                if ss_i2 > ss_i1 && ss_i2 <= numel(disp_um)
                    disp_steady = disp_um(ss_i1:ss_i2);
                    disp_steady = disp_steady - mean(disp_steady);
                    
                    rms_um = sqrt(mean(disp_steady.^2));
                    
                    % True physical peak-to-peak displacement
                    min_dist = max(1, round(0.8 * (Fs / C.FreqHz))); 
                    [pos_peaks, ~] = findpeaks(disp_steady, 'MinPeakDistance', min_dist);
                    [neg_peaks, ~] = findpeaks(-disp_steady, 'MinPeakDistance', min_dist);
                    
                    if ~isempty(pos_peaks) && ~isempty(neg_peaks)
                        pp_um = mean(pos_peaks) - mean(-neg_peaks); 
                    else
                        pp_um = rms_um * 2 * sqrt(2); % Fallback
                    end
                else
                    rms_um = NaN;
                    pp_um  = NaN;
                end

                % --- Store ---
                idx = idx + 1;
                disp_results(idx).Position          = VibData(p).Name;
                disp_results(idx).PosIdx            = p;
                disp_results(idx).PotStep           = C.PotStep;
                disp_results(idx).FreqHz            = C.FreqHz;
                disp_results(idx).BurstIdx          = k;
                disp_results(idx).Displacement_um   = pp_um;
                disp_results(idx).RMS_um            = rms_um;
                disp_results(idx).RMS_g             = rms_g;
                disp_results(idx).SteadyStart_ms    = steady_start_ms;
                disp_results(idx).SteadyDur_ms      = steady_dur_ms;
                disp_results(idx).Source             = 'Calibration';
            end
        end
        fprintf('  %s: %d bursts\n', VibData(p).Name, ...
            sum(strcmp({disp_results.Position}, VibData(p).Name)));
    end

    fprintf('\nTotal: %d displacement measurements.\n', idx);
    if idx == 0, warning('No displacement computed.'); return; end

    %% === SUMMARY TABLES ===
    all_pot  = [disp_results.PotStep];
    all_freq = [disp_results.FreqHz];
    all_disp = [disp_results.Displacement_um];
    all_rms  = [disp_results.RMS_um];
    all_ss   = [disp_results.SteadyStart_ms];

    fprintf('\n=========================================================================\n');
    fprintf('  DISPLACEMENT BY FREQUENCY x POT STEP\n');
    fprintf('=========================================================================\n');
    fprintf('%-8s | %-6s | %-18s | %-18s | %-12s | %-4s\n', ...
            'FREQ', 'POT', 'P-P (um)', 'RMS (um)', 'SS onset (ms)', 'N');
    fprintf('%s\n', repmat('-', 1, 78));

    u_freq = sort(unique(all_freq));
    u_pot  = sort(unique(all_pot));

    for fi = 1:numel(u_freq)
        for pi = 1:numel(u_pot)
            f = u_freq(fi); p = u_pot(pi);
            mask = (all_freq == f) & (all_pot == p);
            if sum(mask) == 0, continue; end
            d = all_disp(mask); r = all_rms(mask); s = all_ss(mask);
            fprintf('%3d Hz   | P%-4d | %6.2f +/- %-7.2f | %6.2f +/- %-7.2f | %5.1f +/- %-4.1f | %d\n', ...
                f, p, mean(d,'omitnan'), std(d,'omitnan'), ...
                      mean(r,'omitnan'), std(r,'omitnan'), ...
                      mean(s,'omitnan'), std(s,'omitnan'), sum(mask));
        end
    end

    %% === DIAGNOSTIC PLOTS (conditional) ===
    if do_plot
        plot_displacement_diagnostics(disp_results, u_freq, u_pot);
    end
end

%% ========================================================================
function plot_displacement_diagnostics(res, u_freq, u_pot)
    all_pot  = [res.PotStep];
    all_freq = [res.FreqHz];
    all_disp = [res.Displacement_um];
    all_pidx = [res.PosIdx];

    % --- Figure 1: Displacement vs Pot Step (one line per freq) ---
    figure('Name', 'Displacement: Transfer Function', 'Color', 'w', ...
           'Position', [100 100 700 500]);
    hold on; grid on; box on;

    colors = lines(numel(u_freq));
    markers = {'o', 's', '^'};

    for fi = 1:numel(u_freq)
        f = u_freq(fi);
        mean_vals = zeros(numel(u_pot), 1);
        std_vals  = zeros(numel(u_pot), 1);

        for pi = 1:numel(u_pot)
            p = u_pot(pi);
            mask = (all_freq == f) & (all_pot == p);
            d = all_disp(mask);
            mean_vals(pi) = mean(d, 'omitnan');
            std_vals(pi)  = std(d, 'omitnan');
        end

        mk = markers{min(fi, numel(markers))};
        errorbar(u_pot, mean_vals, std_vals, ['-' mk], ...
                 'Color', colors(fi,:), 'LineWidth', 1.5, ...
                 'MarkerSize', 8, 'MarkerFaceColor', colors(fi,:), ...
                 'DisplayName', sprintf('%d Hz', f));
    end

    xlabel('Pot Step', 'FontSize', 12);
    ylabel('Peak-to-Peak Displacement (\mum)', 'FontSize', 12);
    legend('Location', 'northwest', 'FontSize', 11);
    set(gca, 'FontSize', 11, 'TickDir', 'out');
    title('Amplitude Transfer Function');

    % --- Figure 2: Per-position consistency ---
    figure('Name', 'Displacement: Position Consistency', 'Color', 'w', ...
           'Position', [820 100 700 500]);
    hold on; grid on; box on;

    u_pos = unique(all_pidx);
    colors_pos = lines(numel(u_pos));

    % Plot at 250 Hz to show spatial uniformity
    target_f = 250;
    for qi = 1:numel(u_pos)
        pos = u_pos(qi);
        mean_vals = zeros(numel(u_pot), 1);

        for pi = 1:numel(u_pot)
            p = u_pot(pi);
            mask = (all_freq == target_f) & (all_pot == p) & (all_pidx == pos);
            d = all_disp(mask);
            mean_vals(pi) = mean(d, 'omitnan');
        end

        plot(u_pot, mean_vals, '-o', 'Color', colors_pos(qi,:), ...
             'LineWidth', 1.2, 'MarkerSize', 6, ...
             'MarkerFaceColor', colors_pos(qi,:), ...
             'DisplayName', sprintf('Pos%d', pos));
    end

    xlabel('Pot Step', 'FontSize', 12);
    ylabel('Peak-to-Peak Displacement (\mum)', 'FontSize', 12);
    legend('Location', 'northwest', 'FontSize', 10);
    set(gca, 'FontSize', 11, 'TickDir', 'out');
    title(sprintf('Spatial Uniformity at %d Hz', target_f));
end