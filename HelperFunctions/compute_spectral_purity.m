function spectral_results = compute_spectral_purity(VibData, do_plot, pos_idx)
% COMPUTE_SPECTRAL_PURITY  Time-frequency analysis of vibrotactile bursts.
%
%   spectral_results = compute_spectral_purity(VibData)
%   spectral_results = compute_spectral_purity(VibData, do_plot)
%   spectral_results = compute_spectral_purity(VibData, do_plot, pos_idx)
%
%   Produces for each condition (single representative trial):
%     Panel A: CWT scalogram  — onset transient zoom (-50 to +200 ms)
%     Panel B: STFT spectrogram — full burst (-100 to +1100 ms for cal,
%              full train for trains)
%     Panel C: Time-resolved THD% curve overlaid on amplitude envelope
%
%   Calibration: one figure per frequency (64, 128, 250 Hz),
%                rows = pot steps, 3 columns = CWT | STFT | THD
%   Trains:      one figure per pot step (3, 18, 35),
%                rows = IBI conditions, 3 columns = CWT | STFT | THD
%
%   Requires: Wavelet Toolbox, Signal Processing Toolbox
%
%   The returned struct (from the quantitative summary) covers calibration
%   bursts at a single position (pos_idx, default 1); the CWT/STFT/THD
%   figures additionally include trains.
%
%   Output struct array fields (steady-state, middle 80% of each burst):
%     .PotStep .FreqHz .Source      — condition identifiers
%     .THD_SteadyState_Pct          — mean THD across steady-state frames (%)
%     .THD_SD_Pct                   — SD of THD across frames (%)
%     .Fund_Power_dB                — mean fundamental power (dB)
%     .Rep                          — repetition index within condition

    if nargin < 2, do_plot = true; end
    if nargin < 3, pos_idx = 1; end

    D = VibData(pos_idx).Raw;
    if isempty(D)
        warning('No data for position %d', pos_idx);
        spectral_results = [];
        return;
    end

    Fs = D.Fs;
    Z  = double(D.Z);

    fprintf('\n================================================================\n');
    fprintf('  SPECTRAL PURITY ANALYSIS — %s (Fs = %.0f Hz)\n', ...
            VibData(pos_idx).Name, Fs);
    fprintf('================================================================\n');

    % =====================================================================
    %  STFT PARAMETERS
    % =====================================================================
    STFT_WIN_MS      = 50;                          % 50 ms window
    STFT_OVERLAP_PCT = 0.90;                        % 90% overlap
    STFT_NFFT        = 1024;                        % zero-pad for smooth freq axis
    FREQ_MAX         = 800;                         % Hz display limit

    stft_win_samp    = round(STFT_WIN_MS / 1000 * Fs);
    stft_overlap     = round(stft_win_samp * STFT_OVERLAP_PCT);
    stft_window      = hann(stft_win_samp);

    % CWT onset zoom window
    CWT_PRE_MS  = 50;
    CWT_POST_MS = 200;

    % STFT full-burst window
    STFT_PRE_MS      = 100;
    STFT_POST_CAL_MS = 1100;   % calibration: 1000ms burst + 100ms tail

    % THD settings
    THD_HARM_COUNT   = 6;      % track up to 6th harmonic
    THD_BIN_TOL_HZ   = 8;      % +/- Hz around each harmonic peak

    % =====================================================================
    %  SECTION 1: CALIBRATION BURSTS (figures only if do_plot)
    % =====================================================================
    if do_plot && isfield(VibData(pos_idx), 'Calibration') && ~isempty(VibData(pos_idx).Calibration)
        Cal = VibData(pos_idx).Calibration;

        % Group by frequency
        cal_freqs = unique([Cal.FreqHz]);
        cal_freqs = sort(cal_freqs);  % 64, 128, 250

        for fi = 1:numel(cal_freqs)
            freq = cal_freqs(fi);

            % Find all conditions at this frequency
            freq_mask = [Cal.FreqHz] == freq;
            conds = Cal(freq_mask);
            n_conds = numel(conds);

            % One figure: rows = pot steps, cols = CWT | STFT | THD
            figure('Color', 'w', ...
                     'Name', sprintf('Spectral: Cal %d Hz (%s)', freq, VibData(pos_idx).Name), ...
                     'Position', [40+fi*30, 40+fi*30, 1500, 220*n_conds + 60]);

            for ci = 1:n_conds
                C = conds(ci);
                if C.NBursts == 0, continue; end

                % Take first rep as representative
                onset = C.OnsetIdx(1);
                pot   = C.PotStep;

                % --- Extract segments ---
                [seg_cwt, t_cwt]   = extract_segment(Z, onset, Fs, CWT_PRE_MS, CWT_POST_MS);
                [seg_stft, t_stft] = extract_segment(Z, onset, Fs, STFT_PRE_MS, STFT_POST_CAL_MS);

                if isempty(seg_cwt) || isempty(seg_stft), continue; end

                % --- Panel A: CWT scalogram (onset zoom) ---
                ax_cwt = subplot(n_conds, 3, (ci-1)*3 + 1);
                plot_cwt_panel(ax_cwt, seg_cwt, Fs, t_cwt, FREQ_MAX, ...
                               sprintf('CWT: Pot %d, %d Hz', pot, freq));

                % --- Panel B: STFT spectrogram (full burst) ---
                ax_stft = subplot(n_conds, 3, (ci-1)*3 + 2);
                plot_stft_panel(ax_stft, seg_stft, Fs, t_stft, ...
                                stft_window, stft_overlap, STFT_NFFT, FREQ_MAX, ...
                                sprintf('STFT: Pot %d, %d Hz', pot, freq));

                % --- Panel C: THD curve ---
                ax_thd = subplot(n_conds, 3, (ci-1)*3 + 3);
                plot_thd_panel(ax_thd, seg_stft, Fs, t_stft, freq, ...
                               stft_window, stft_overlap, STFT_NFFT, ...
                               THD_HARM_COUNT, THD_BIN_TOL_HZ, ...
                               sprintf('THD: Pot %d, %d Hz', pot, freq));
            end

            sgtitle(sprintf('Spectral Characterisation: %d Hz — %s (rep 1)', ...
                    freq, VibData(pos_idx).Name), ...
                    'FontSize', 14, 'FontWeight', 'bold');
        end
        fprintf('  Calibration spectral plots done.\n');
    end

    % =====================================================================
    %  SECTION 2: TRAINS (figures only if do_plot)
    % =====================================================================
    if do_plot && isfield(VibData(pos_idx), 'Trains') && ~isempty(VibData(pos_idx).Trains)
        Tr = VibData(pos_idx).Trains;

        train_pots = unique([Tr.PotStep]);
        train_ibis = unique([Tr.IBI_ms]);

        for pi = 1:numel(train_pots)
            pot = train_pots(pi);

            % Find IBI conditions at this pot (take rep 1)
            pot_mask = ([Tr.PotStep] == pot) & ([Tr.RepIdx] == 1);
            pot_trains = Tr(pot_mask);

            if isempty(pot_trains), continue; end

            % Sort by IBI
            [~, si] = sort([pot_trains.IBI_ms]);
            pot_trains = pot_trains(si);
            n_ibi = numel(pot_trains);

            fig = figure('Color', 'w', ...
                         'Name', sprintf('Spectral: Train Pot %d (%s)', pot, VibData(pos_idx).Name), ...
                         'Position', [60+pi*30, 60+pi*30, 1500, 220*n_ibi + 60]);

            for ii = 1:n_ibi
                T = pot_trains(ii);
                if T.NBursts == 0, continue; end

                % Calculate train duration
                first_onset = T.BurstOnsetIdx(1);
                last_onset  = T.BurstOnsetIdx(end);
                burst_samps = round(T.BurstMs / 1000 * Fs);
                train_end   = last_onset + burst_samps;
                train_dur_ms = (train_end - first_onset) / Fs * 1000;

                % Windows
                train_stft_post = train_dur_ms + 400;  % 400ms tail after last burst
                train_cwt_post  = min(200, train_dur_ms);  % zoom on first burst(s)

                % --- Extract segments ---
                [seg_cwt, t_cwt]   = extract_segment(Z, first_onset, Fs, CWT_PRE_MS, train_cwt_post);
                [seg_stft, t_stft] = extract_segment(Z, first_onset, Fs, STFT_PRE_MS, train_stft_post);

                if isempty(seg_cwt) || isempty(seg_stft), continue; end

                ibi_str = sprintf('%.1f', T.IBI_ms);

                % --- Panel A: CWT ---
                ax_cwt = subplot(n_ibi, 3, (ii-1)*3 + 1);
                plot_cwt_panel(ax_cwt, seg_cwt, Fs, t_cwt, FREQ_MAX, ...
                               sprintf('CWT: P%d IBI %s ms', pot, ibi_str));

                % --- Panel B: STFT ---
                ax_stft = subplot(n_ibi, 3, (ii-1)*3 + 2);
                plot_stft_panel(ax_stft, seg_stft, Fs, t_stft, ...
                                stft_window, stft_overlap, STFT_NFFT, FREQ_MAX, ...
                                sprintf('STFT: P%d IBI %s ms', pot, ibi_str));

                % --- Panel C: THD ---
                ax_thd = subplot(n_ibi, 3, (ii-1)*3 + 3);
                plot_thd_panel(ax_thd, seg_stft, Fs, t_stft, T.FreqHz, ...
                               stft_window, stft_overlap, STFT_NFFT, ...
                               THD_HARM_COUNT, THD_BIN_TOL_HZ, ...
                               sprintf('THD: P%d IBI %s ms', pot, ibi_str));
            end

            sgtitle(sprintf('Train Spectral: Pot %d, 250 Hz — %s (rep 1)', ...
                    pot, VibData(pos_idx).Name), ...
                    'FontSize', 14, 'FontWeight', 'bold');
        end
        fprintf('  Train spectral plots done.\n');
    end

    % =====================================================================
    %  QUANTITATIVE SUMMARY
    % =====================================================================
    spectral_results = compute_thd_summary(VibData, pos_idx, ...
                                           stft_window, stft_overlap, STFT_NFFT, ...
                                           THD_HARM_COUNT, THD_BIN_TOL_HZ);
end

%% ========================================================================
%  SEGMENT EXTRACTION
%  ========================================================================
function [seg, t_ms] = extract_segment(Z, onset_idx, Fs, pre_ms, post_ms)
    pre_samp  = round(pre_ms / 1000 * Fs);
    post_samp = round(post_ms / 1000 * Fs);

    i1 = onset_idx - pre_samp;
    i2 = onset_idx + post_samp;

    if i1 < 1 || i2 > numel(Z)
        seg = []; t_ms = [];
        return;
    end

    seg  = Z(i1:i2);
    seg  = seg - mean(seg);  % remove DC
    t_ms = ((0:numel(seg)-1) - pre_samp) / Fs * 1000;
end

%% ========================================================================
%  PANEL A: CWT SCALOGRAM
%  ========================================================================
function plot_cwt_panel(ax, seg, Fs, t_ms, freq_max, title_str)
    axes(ax); %#ok<LAXES>

    % CWT with analytic Morse wavelet
    [wt, f_cwt] = cwt(seg, Fs);

    % Convert to power (dB)
    P = 10 * log10(abs(wt).^2 + eps);

    % Frequency mask
    f_mask = f_cwt <= freq_max;
    P_plot = P(f_mask, :);
    f_plot = f_cwt(f_mask);

    % Note: cwt returns frequencies high-to-low, so flip for correct orientation
    imagesc(t_ms, f_plot, P_plot);
    axis xy;  % frequency ascending

    colormap(ax, parula);
    cb = colorbar;
    cb.Label.String = 'dB';
    cb.FontSize = 7;

    % TTL marker
    hold on;
    xline(0, 'w--', 'LineWidth', 1.5);

    xlabel('Time (ms)');
    ylabel('Freq (Hz)');
    title(title_str, 'FontSize', 9, 'FontWeight', 'bold');
    ylim([10 freq_max]);
    set(ax, 'FontSize', 8);
end

%% ========================================================================
%  PANEL B: STFT SPECTROGRAM
%  ========================================================================
function plot_stft_panel(ax, seg, Fs, t_ms, win, noverlap, nfft, freq_max, title_str)
    axes(ax); %#ok<LAXES>

    [S, f_stft, t_stft] = spectrogram(seg, win, noverlap, nfft, Fs);

    % Convert spectrogram time to ms relative to onset
    % t_stft is in seconds from start of segment
    t_stft_ms = t_stft * 1000 + t_ms(1);  % offset by segment pre-time

    P = 10 * log10(abs(S).^2 + eps);

    f_mask = f_stft <= freq_max;

    imagesc(t_stft_ms, f_stft(f_mask), P(f_mask, :));
    axis xy;
    colormap(ax, parula);
    cb = colorbar;
    cb.Label.String = 'dB';
    cb.FontSize = 7;

    hold on;
    xline(0, 'w--', 'LineWidth', 1.5);

    xlabel('Time (ms)');
    ylabel('Freq (Hz)');
    title(title_str, 'FontSize', 9, 'FontWeight', 'bold');
    ylim([0 freq_max]);
    set(ax, 'FontSize', 8);
end

%% ========================================================================
%  PANEL C: TIME-RESOLVED THD
%  ========================================================================
function plot_thd_panel(ax, seg, Fs, t_ms, fund_freq, win, noverlap, nfft, ...
                        n_harmonics, bin_tol_hz, title_str)
    axes(ax); %#ok<LAXES>

    [S, f_stft, t_stft] = spectrogram(seg, win, noverlap, nfft, Fs);
    t_stft_ms = t_stft * 1000 + t_ms(1);

    mag = abs(S);  % magnitude spectrum at each time frame
    df  = f_stft(2) - f_stft(1);  % frequency resolution

    n_frames = size(mag, 2);
    thd_pct       = NaN(1, n_frames);
    fund_power_db = NaN(1, n_frames);

    for k = 1:n_frames
        spectrum = mag(:, k);

        % Fundamental power: sum in band around f0
        f0_band = abs(f_stft - fund_freq) <= bin_tol_hz;
        P_fund  = sum(spectrum(f0_band).^2);

        % Harmonic power: sum across harmonics 2..N
        P_harm = 0;
        for h = 2:(n_harmonics + 1)
            fh = fund_freq * h;
            if fh > Fs/2, break; end
            h_band = abs(f_stft - fh) <= bin_tol_hz;
            P_harm = P_harm + sum(spectrum(h_band).^2);
        end

        if P_fund > 0
            thd_pct(k) = sqrt(P_harm / P_fund) * 100;
        end
        fund_power_db(k) = 10 * log10(P_fund + eps);
    end

    % --- Plot ---
    yyaxis(ax, 'left');
    plot(t_stft_ms, thd_pct, 'r-', 'LineWidth', 1.2);
    ylabel('THD (%)');
    ylim_upper = min(200, prctile(thd_pct(~isnan(thd_pct)), 98) * 1.3);
    if isnan(ylim_upper) || ylim_upper <= 0, ylim_upper = 100; end
    ylim([0 ylim_upper]);
    ax.YColor = 'r';

    yyaxis(ax, 'right');
    plot(t_stft_ms, fund_power_db, 'b-', 'LineWidth', 1.0);
    ylabel(sprintf('Fund. Power (dB) @ %d Hz', fund_freq));
    ax.YColor = 'b';

    hold on;
    xline(0, 'k--', 'LineWidth', 1.5);

    xlabel('Time (ms)');
    title(title_str, 'FontSize', 9, 'FontWeight', 'bold');
    set(ax, 'FontSize', 8);

    % Reset to left axis as default
    yyaxis(ax, 'left');
end

%% ========================================================================
%  QUANTITATIVE THD SUMMARY (all reps, steady-state only)
%  ========================================================================
function results = compute_thd_summary(VibData, pos_idx, win, noverlap, nfft, ...
                                       n_harmonics, bin_tol_hz)
    D = VibData(pos_idx).Raw;
    Fs = D.Fs;
    Z  = double(D.Z);

    results = struct('PotStep',{}, 'FreqHz',{}, 'Source',{}, ...
                     'THD_SteadyState_Pct',{}, 'THD_SD_Pct',{}, ...
                     'Fund_Power_dB',{}, 'Rep',{});
    idx = 0;

    % --- Calibration ---
    if isfield(VibData(pos_idx), 'Calibration') && ~isempty(VibData(pos_idx).Calibration)
        Cal = VibData(pos_idx).Calibration;
        for ci = 1:numel(Cal)
            C = Cal(ci);
            for k = 1:C.NBursts
                onset = C.OnsetIdx(k);
                burst_ms = C.BurstMs;

                % Extract steady-state middle 80%
                margin_ms = burst_ms * 0.1;
                ss_start  = onset + round(margin_ms / 1000 * Fs);
                ss_end    = onset + round((burst_ms - margin_ms) / 1000 * Fs);
                if ss_end > numel(Z), continue; end

                seg = Z(ss_start:ss_end) - mean(Z(ss_start:ss_end));

                [S, f_stft] = spectrogram(seg, win, noverlap, nfft, Fs);
                mag = abs(S);

                % Average THD across all frames in steady state
                thd_frames = compute_thd_frames(mag, f_stft, C.FreqHz, Fs, ...
                                                n_harmonics, bin_tol_hz);

                idx = idx + 1;
                results(idx).PotStep            = C.PotStep;
                results(idx).FreqHz             = C.FreqHz;
                results(idx).Source             = 'Calibration';
                results(idx).THD_SteadyState_Pct = mean(thd_frames, 'omitnan');
                results(idx).THD_SD_Pct         = std(thd_frames, 'omitnan');
                results(idx).Fund_Power_dB      = mean(10*log10(compute_fund_power(mag, f_stft, C.FreqHz, bin_tol_hz) + eps));
                results(idx).Rep                = k;
            end
        end
    end

    % --- Print summary ---
    if idx == 0
        fprintf('  No spectral results computed.\n');
        return;
    end

    all_pot  = [results.PotStep];
    all_freq = [results.FreqHz];
    all_thd  = [results.THD_SteadyState_Pct];
    all_fund = [results.Fund_Power_dB];

    u_freq = sort(unique(all_freq));
    u_pot  = unique(all_pot);

    fprintf('\n=========================================================================\n');
    fprintf('  STEADY-STATE THD SUMMARY (%s, middle 80%%)\n', VibData(pos_idx).Name);
    fprintf('=========================================================================\n');
    fprintf('%-8s | %-6s | %-15s | %-15s | %-6s\n', ...
            'Freq', 'Pot', 'THD (%)', 'Fund Power(dB)', 'N');
    fprintf('%s\n', repmat('-', 1, 62));

    for fi = 1:numel(u_freq)
        freq = u_freq(fi);
        for pi = 1:numel(u_pot)
            pot = u_pot(pi);
            mask = (all_freq == freq) & (all_pot == pot);
            if ~any(mask), continue; end
            fprintf('%3d Hz   | P%-4d | %5.1f +/- %-6.1f | %5.1f +/- %-6.1f | %d\n', ...
                freq, pot, ...
                mean(all_thd(mask)), std(all_thd(mask)), ...
                mean(all_fund(mask)), std(all_fund(mask)), ...
                sum(mask));
        end
    end
    fprintf('=========================================================================\n\n');
end

%% ========================================================================
%  THD helper: compute THD for each frame of a spectrogram magnitude matrix
%  ========================================================================
function thd_frames = compute_thd_frames(mag, f_stft, fund_freq, Fs, n_harmonics, bin_tol)
    n_frames = size(mag, 2);
    thd_frames = NaN(1, n_frames);

    for k = 1:n_frames
        spectrum = mag(:, k);

        f0_band = abs(f_stft - fund_freq) <= bin_tol;
        P_fund  = sum(spectrum(f0_band).^2);

        P_harm = 0;
        for h = 2:(n_harmonics + 1)
            fh = fund_freq * h;
            if fh > Fs/2, break; end
            h_band = abs(f_stft - fh) <= bin_tol;
            P_harm = P_harm + sum(spectrum(h_band).^2);
        end

        if P_fund > 0
            thd_frames(k) = sqrt(P_harm / P_fund) * 100;
        end
    end
end

%% ========================================================================
%  Fund power helper
%  ========================================================================
function P_fund = compute_fund_power(mag, f_stft, fund_freq, bin_tol)
    n_frames = size(mag, 2);
    P_fund = zeros(1, n_frames);
    for k = 1:n_frames
        f0_band = abs(f_stft - fund_freq) <= bin_tol;
        P_fund(k) = sum(mag(f0_band, k).^2);
    end
end