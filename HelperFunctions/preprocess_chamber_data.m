function VibData = preprocess_chamber_data(Folder)
% PREPROCESS_CHAMBER_DATA  Import & segment unified v3 chamber recordings.
%
%   VibData = preprocess_chamber_data(Folder)
%
%   Protocol (Arduino v3 — unified calibration):
%
%   PREAMBLE:     5 x 100ms TTL marker pulses, then 1s baseline
%   CALIBRATION:  pot [1,3,5,10,18,35,61] x freq [64,128,250] x 5 reps
%                 1000ms ON, 4s ITI between all bursts
%                 4s pot-change gap with TTL pot-marker packet (mid-gap pot step)
%                 Order: pot(outer) -> freq(middle) -> rep(inner)
%   TRAINS:       pot [3,18,35] x IBI [27.5,33.35,66.67,167] x 5 reps
%                 6 x 100ms bursts/train, 4s ITI between trains
%                 4s pot-change gap between pot groups
%                 Order: pot(outer) -> IBI(middle) -> rep(inner)
%
%   TTL pot-marker packets (100ms header + N x 30ms data pulses) fire
%   at the start of each pot-change gap.  These are automatically
%   filtered from the burst stream.
%
%   Output: VibData(1..5) struct array with fields:
%     .Name                -- 'Pos1'..'Pos5'
%     .Raw                 -- struct with t, X, Y, Z, TTL, Fs, seq
%     .Markers             -- onset indices of preamble marker pulses
%     .Calibration(1..21)  -- 7 pot x 3 freq conditions (5 reps each)
%     .Trains(1..60)       -- 3 pot x 4 IBI x 5 reps

%% === GROUND TRUTH PARAMETERS ===
filePrefixes         = {'Pos1', 'Pos2', 'Pos3', 'Pos4', 'Pos5'};
filePrefixesFallback = {'OST',  'SYD',  'VEST', 'NORD', 'MID'};
nPos = 5;

% --- Calibration (unified, all freqs at every pot) ---
CAL_FREQS      = [64, 128, 250];       % Hz — order within each pot block
CAL_N_FREQ     = 3;
CAL_POT_STEPS  = [1, 3, 5, 10, 18, 35, 61];
CAL_N_POT      = 7;
CAL_REPS       = 5;                    % per (pot, freq) pair
CAL_BURST_MS   = 1000;                 % 1 s ON
CAL_TOTAL      = CAL_N_POT * CAL_N_FREQ * CAL_REPS;   % 105

% --- Trains ---
TRAIN_FREQ       = 250;                % Hz (fixed)
TRAIN_POT_STEPS  = [3, 18, 35];
TRAIN_N_POT      = 3;
TRAIN_IBI_MS     = [27.5, 33.35, 66.67, 167];   % 4 IBI conditions
TRAIN_N_IBI      = 4;
TRAIN_REPS       = 5;
TRAIN_BURSTS_PER = 6;                  % bursts per train
TRAIN_BURST_MS   = 100;                % ms per burst
TRAIN_TOTAL_TRAINS  = TRAIN_N_POT * TRAIN_N_IBI * TRAIN_REPS;   % 60
TRAIN_TOTAL_BURSTS  = TRAIN_TOTAL_TRAINS * TRAIN_BURSTS_PER;    % 360

% --- Detection thresholds ---
CAL_DUR_THRESH_MS         = 500;    % >500ms = calibration burst
POT_MARKER_BIT_MAX_MS     = 50;     % pot-marker data pulses ~30ms
POT_MARKER_BIT_MIN_MS     = 10;     % noise floor
TRAIN_BOUNDARY_THRESH_MS  = 800;    % max IBI (167) < 800 < ITI (4000)
BASELINE_GAP_THRESH_S     = 1.5;    % gap after preamble markers
FFT_FREQ_TOLERANCE        = 15;     % Hz
MIN_CLUSTER_BURSTS        = 2;      % clusters with <2 bursts = pot-marker headers

%% === 1. LOAD RAW DATA ===
VibData = struct();
for i = 1:nPos
    VibData(i).Name = sprintf('Pos%d', i);
    fpath = fullfile(Folder, [filePrefixes{i} '_Accelerometers.bin']);
    if ~isfile(fpath)
        fpath = fullfile(Folder, [filePrefixesFallback{i} '_Accelerometers.bin']);
    end
    if ~isfile(fpath)
        warning('File not found: %s -- skipping %s', fpath, VibData(i).Name);
        VibData(i).Raw = [];
        continue;
    end
    fprintf('Loading %s ...\n', VibData(i).Name);
    VibData(i).Raw = read_adxl_bin(fpath);
end
fprintf('Import complete.\n\n');

%% === 2. SEGMENT EACH POSITION ===
for i = 1:nPos
    D = VibData(i).Raw;
    if isempty(D)
        VibData(i).Calibration = [];
        VibData(i).Trains      = [];
        VibData(i).Markers     = [];
        continue;
    end

    fprintf('===== %s =====\n', VibData(i).Name);

    % --- Detect all edges ---
    [all_onsets, all_times] = detect_all_rising_edges(D.t, D.TTL);
    all_falls = find(diff(double(D.TTL)) == -1);
    fprintf('  Total rising edges: %d\n', numel(all_onsets));

    % --- Measure burst durations ---
    burst_dur_ms = measure_burst_durations(D.t, all_onsets, all_falls);

    % --- Separate preamble markers from real content ---
    first_real_idx = 1;
    marker_count   = 0;
    for k = 2:numel(all_onsets)
        prev_fall = max(all_falls(all_falls < all_onsets(k)));
        if isempty(prev_fall), continue; end
        gap = D.t(all_onsets(k)) - D.t(prev_fall);
        if gap > BASELINE_GAP_THRESH_S
            first_real_idx = k;
            marker_count   = k - 1;
            break;
        end
    end

    fprintf('  Preamble markers: %d', marker_count);
    if marker_count > 0
        fprintf(' (durations:');
        fprintf(' %.0f', burst_dur_ms(1:marker_count));
        fprintf(' ms)');
    end
    fprintf('\n');
    VibData(i).Markers = all_onsets(1:marker_count);

    % --- Work with post-preamble edges only ---
    real_onsets  = all_onsets(first_real_idx:end);
    real_times   = all_times(first_real_idx:end);
    real_dur_ms  = burst_dur_ms(first_real_idx:end);

    % =================================================================
    %  STEP A: Remove pot-marker DATA pulses (~30ms)
    %  These are the short data bits in the TTL marker packets.
    %  The 100ms header pulses are handled later by cluster filtering.
    % =================================================================
    is_marker_bit = (real_dur_ms >= POT_MARKER_BIT_MIN_MS) & ...
                    (real_dur_ms <  POT_MARKER_BIT_MAX_MS);
    n_marker_bits = sum(is_marker_bit);
    fprintf('  Pot-marker data pulses removed: %d\n', n_marker_bits);

    filt_onsets  = real_onsets(~is_marker_bit);
    filt_times   = real_times(~is_marker_bit);
    filt_dur_ms  = real_dur_ms(~is_marker_bit);

    % =================================================================
    %  STEP B: Classify by duration
    %    > 500ms  → calibration burst  (1000ms nominal)
    %    ≤ 500ms  → candidate train burst OR pot-marker header (100ms)
    % =================================================================
    is_cal = filt_dur_ms > CAL_DUR_THRESH_MS;

    cal_onsets  = filt_onsets(is_cal);
    cal_times   = filt_times(is_cal);

    short_onsets  = filt_onsets(~is_cal);
    short_times   = filt_times(~is_cal);
    short_dur_ms  = filt_dur_ms(~is_cal);

    n_cal   = numel(cal_onsets);
    n_short = numel(short_onsets);
    fprintf('  Calibration bursts (>500ms): %d (expected %d)\n', n_cal, CAL_TOTAL);

    % =================================================================
    %  STEP C: Cluster short pulses into trains
    %  Gap > 800ms marks a new cluster.
    %  Clusters with <2 bursts = isolated pot-marker headers → discard.
    % =================================================================
    if n_short > 0
        short_gaps_ms = NaN(n_short, 1);
        for k = 2:n_short
            short_gaps_ms(k) = (short_times(k) - short_times(k-1)) * 1000;
        end

        is_cluster_start = true(n_short, 1);
        if n_short > 1
            is_cluster_start(2:end) = short_gaps_ms(2:end) > TRAIN_BOUNDARY_THRESH_MS;
        end

        cluster_starts = find(is_cluster_start);
        n_clusters     = numel(cluster_starts);
        cluster_ends   = [cluster_starts(2:end) - 1; n_short];
        cluster_sizes  = cluster_ends - cluster_starts + 1;

        % Filter: keep only real trains (≥2 bursts)
        real_mask = cluster_sizes >= MIN_CLUSTER_BURSTS;
        n_pot_headers = sum(~real_mask);
        fprintf('  Pot-marker headers removed: %d (isolated short pulses)\n', n_pot_headers);

        train_cl_starts = cluster_starts(real_mask);
        train_cl_ends   = cluster_ends(real_mask);
        train_cl_sizes  = cluster_sizes(real_mask);
        n_trains_detected = numel(train_cl_starts);

        % Flatten real train bursts for counting
        real_train_burst_mask = false(n_short, 1);
        for ci = 1:n_trains_detected
            real_train_burst_mask(train_cl_starts(ci):train_cl_ends(ci)) = true;
        end
        train_onsets = short_onsets(real_train_burst_mask);
        train_times  = short_times(real_train_burst_mask);
        n_train_bursts = numel(train_onsets);
    else
        n_trains_detected = 0;
        train_onsets = [];
        train_times  = [];
        n_train_bursts = 0;
        train_cl_starts = [];
        train_cl_ends   = [];
        train_cl_sizes  = [];
    end

    fprintf('  Train clusters: %d (expected %d)\n', n_trains_detected, TRAIN_TOTAL_TRAINS);
    fprintf('  Train bursts:   %d (expected %d)\n', n_train_bursts, TRAIN_TOTAL_BURSTS);

    % =================================================================
    %  CALIBRATION: FFT each burst → frequency, then assign sequentially
    %  Order: pot(outer) -> freq(middle) -> rep(inner)
    %  Per pot block: 64 x5, 128 x5, 250 x5 = 15 bursts
    % =================================================================
    cal_freqs = zeros(n_cal, 1);
    for k = 1:n_cal
        cal_freqs(k) = estimate_burst_frequency(D, cal_onsets(k), CAL_BURST_MS);
    end

    Calibration = struct();
    cond_idx    = 0;
    burst_cursor = 0;

    for p = 1:CAL_N_POT
        for f = 1:CAL_N_FREQ
            cond_idx  = cond_idx + 1;
            idx_range = (1:CAL_REPS) + burst_cursor;
            idx_range = idx_range(idx_range <= n_cal);

            Calibration(cond_idx).PotStep       = CAL_POT_STEPS(p);
            Calibration(cond_idx).FreqHz        = CAL_FREQS(f);
            Calibration(cond_idx).BurstMs       = CAL_BURST_MS;
            Calibration(cond_idx).OnsetIdx      = cal_onsets(idx_range);
            Calibration(cond_idx).OnsetTimes    = cal_times(idx_range);
            Calibration(cond_idx).NBursts       = numel(idx_range);

            if ~isempty(idx_range)
                Calibration(cond_idx).MeasuredFreqs = cal_freqs(idx_range);
            else
                Calibration(cond_idx).MeasuredFreqs = [];
            end

            burst_cursor = burst_cursor + CAL_REPS;
        end
    end
    VibData(i).Calibration = Calibration;

    % =================================================================
    %  TRAINS: assign clusters sequentially
    %  Order: pot(outer) -> IBI(middle) -> rep(inner)
    % =================================================================
    Trains = struct();
    train_counter = 0;

    for p = 1:TRAIN_N_POT
        for b = 1:TRAIN_N_IBI
            for r = 1:TRAIN_REPS
                train_counter = train_counter + 1;

                Trains(train_counter).PotStep  = TRAIN_POT_STEPS(p);
                Trains(train_counter).IBI_ms   = TRAIN_IBI_MS(b);
                Trains(train_counter).RepIdx   = r;
                Trains(train_counter).FreqHz   = TRAIN_FREQ;
                Trains(train_counter).BurstMs  = TRAIN_BURST_MS;

                if train_counter <= n_trains_detected
                    ci = train_counter;
                    burst_range_idx = train_cl_starts(ci):train_cl_ends(ci);
                    % Limit to expected bursts per train
                    burst_range_idx = burst_range_idx(1:min(end, TRAIN_BURSTS_PER));

                    Trains(train_counter).TrainOnsetIdx    = short_onsets(burst_range_idx(1));
                    Trains(train_counter).TrainOnsetTime   = short_times(burst_range_idx(1));
                    Trains(train_counter).BurstOnsetIdx    = short_onsets(burst_range_idx);
                    Trains(train_counter).BurstOnsetTimes  = short_times(burst_range_idx);
                    Trains(train_counter).NBursts          = numel(burst_range_idx);
                else
                    Trains(train_counter).TrainOnsetIdx    = [];
                    Trains(train_counter).TrainOnsetTime   = [];
                    Trains(train_counter).BurstOnsetIdx    = [];
                    Trains(train_counter).BurstOnsetTimes  = [];
                    Trains(train_counter).NBursts          = 0;
                end
            end
        end
    end
    VibData(i).Trains = Trains;

    % =================================================================
    %  VALIDATION PRINTOUT
    % =================================================================
    fprintf('\n  --- Validation ---\n');

    fprintf('  CALIBRATION (pot x freq x 5 reps):\n');
    for c = 1:numel(Calibration)
        C = Calibration(c);
        if C.NBursts > 0
            fft_mean = mean(C.MeasuredFreqs, 'omitnan');
            fft_ok   = abs(fft_mean - C.FreqHz) < FFT_FREQ_TOLERANCE;
            fprintf('    Pot %2d @ %3d Hz: %d bursts, FFT = %.0f Hz %s\n', ...
                C.PotStep, C.FreqHz, C.NBursts, fft_mean, ...
                conditional_str(fft_ok, 'OK', '*** MISMATCH ***'));
        else
            fprintf('    Pot %2d @ %3d Hz: 0 bursts *** MISSING ***\n', ...
                C.PotStep, C.FreqHz);
        end
    end

    fprintf('  TRAINS (%d pot x %d IBI x %d reps):\n', ...
        TRAIN_N_POT, TRAIN_N_IBI, TRAIN_REPS);
    n_ok   = sum([Trains.NBursts] == TRAIN_BURSTS_PER);
    n_part = sum([Trains.NBursts] > 0 & [Trains.NBursts] < TRAIN_BURSTS_PER);
    n_miss = sum([Trains.NBursts] == 0);
    fprintf('    %d complete (%d bursts), %d partial, %d missing\n', ...
        n_ok, TRAIN_BURSTS_PER, n_part, n_miss);
    fprintf('\n');
end

%% === 3. SUMMARY TABLE ===
fprintf('================ PROTOCOL SUMMARY (v3) ================\n');
fprintf('%-6s  Cal(%d)  Trains(%d)  TrainBursts(%d)\n', ...
    'Pos', CAL_TOTAL, TRAIN_TOTAL_TRAINS, TRAIN_TOTAL_BURSTS);
for i = 1:nPos
    if isempty(VibData(i).Raw), continue; end
    n_cal = sum([VibData(i).Calibration.NBursts]);
    n_tr  = sum([VibData(i).Trains.NBursts] > 0);
    n_tb  = sum([VibData(i).Trains.NBursts]);
    flag  = '';
    if n_cal ~= CAL_TOTAL || n_tr ~= TRAIN_TOTAL_TRAINS
        flag = ' ***';
    end
    fprintf('%-6s  %4d    %4d        %4d%s\n', ...
        VibData(i).Name, n_cal, n_tr, n_tb, flag);
end
fprintf('========================================================\n\n');
end

%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

function s = conditional_str(cond, s_true, s_false)
    if cond, s = s_true; else, s = s_false; end
end

function [onsets, onset_times] = detect_all_rising_edges(t, TTL)
% DETECT_ALL_RISING_EDGES  Find sample indices where TTL goes 0->1.
    edges  = diff(double(TTL)) == 1;
    onsets = find(edges) + 1;          % +1: edge is between (k) and (k+1)
    onset_times = t(onsets);
end

function burst_dur_ms = measure_burst_durations(t, onsets, falls)
    burst_dur_ms = NaN(numel(onsets), 1);
    for k = 1:numel(onsets)
        next_fall = min(falls(falls > onsets(k)));
        if ~isempty(next_fall)
            burst_dur_ms(k) = (t(next_fall) - t(onsets(k))) * 1000;
        end
    end
end

function peak_freq = estimate_burst_frequency(data, onset_idx, burst_ms)
    Fs = data.Fs;
    burst_samp = round(burst_ms / 1000 * Fs);
    margin = round(0.1 * burst_samp);       % skip first/last 10%
    i1 = onset_idx + margin;
    i2 = onset_idx + burst_samp - margin;
    i2 = min(i2, numel(data.Z));
    if i1 >= i2 || i1 < 1
        peak_freq = NaN;
        return;
    end
    seg = data.Z(i1:i2);
    seg = seg - mean(seg);
    N = length(seg);
    f = (0:N-1) * Fs / N;
    Y = abs(fft(seg));
    Y = Y(1:floor(N/2)+1);
    f = f(1:floor(N/2)+1);
    band = (f >= 20) & (f <= 500);
    [~, idx] = max(Y(band));
    f_band = f(band);
    peak_freq = f_band(idx);
end

function data = read_adxl_bin(filepath)
% READ_ADXL_BIN  Parse ADXL357 binary file (23-byte frames) into struct.
%
%   data = read_adxl_bin(filepath)
%
%   Frame format: [A5 5A] [seq(4)] [t_us(4)] [XYZ(9)] [flags(1)] [CRC(2)] [55]
%                  1  2     3-6       7-10      11-19     20         21-22    23
%
%   Output struct fields:
%     .t    -- time vector (seconds, zeroed to first sample)
%     .X    -- X acceleration (g)
%     .Y    -- Y acceleration (g)
%     .Z    -- Z acceleration (g)
%     .TTL  -- logical: stimulus ON flag (bit 0 of flags byte)
%     .Fs   -- sampling frequency (Hz, from median dt)
%     .seq  -- sequence numbers (uint32)
%     .nFrames -- total valid frames

    data = [];

    if ~isfile(filepath)
        warning('File not found: %s', filepath);
        return;
    end

    fid = fopen(filepath, 'rb');
    raw = fread(fid, Inf, '*uint8');
    fclose(fid);

    FRAME_LEN = 23;
    nBytes = numel(raw);

    if nBytes < FRAME_LEN
        warning('File too small: %s', filepath);
        return;
    end

    % Vectorized frame search
    h1 = raw(1:end-22) == 0xA5;
    h2 = raw(2:end-21) == 0x5A;
    tl = raw(23:end)    == 0x55;
    cand = find(h1 & h2);
    cand = cand(cand + 22 <= nBytes);
    cand = cand(tl(cand));

    if isempty(cand)
        warning('No valid frames found in: %s', filepath);
        return;
    end

    idx_mat = cand(:) + (0:22);
    frames  = reshape(raw(idx_mat.'), FRAME_LEN, []).';
    fprintf('  %s: %d candidate frames', filepath, size(frames,1));

    % CRC16 validation
    CRC_TABLE = build_crc_table();
    n = size(frames, 1);
    crc_calc = repmat(uint16(65535), n, 1);

    for col = 3:20
        byte    = uint16(frames(:, col));
        lut_idx = bitxor(bitand(crc_calc, uint16(255)), byte) + 1;
        crc_calc = bitxor(bitshift(crc_calc, -8), CRC_TABLE(lut_idx));
    end

    crc_file = bitor(bitshift(uint16(frames(:,21)), 8), uint16(frames(:,22)));
    valid    = (crc_calc == crc_file);
    frames   = frames(valid, :);
    fprintf(' -> %d valid (%.1f%%)\n', size(frames,1), 100*sum(valid)/n);

    if isempty(frames)
        warning('No CRC-valid frames in: %s', filepath);
        return;
    end

    % Timestamp (bytes 7-10, big-endian uint32 us)
    t_us   = double(frames(:,7)) * 2^24 + double(frames(:,8)) * 2^16 + ...
             double(frames(:,9)) * 2^8  + double(frames(:,10));
    data.t = (t_us - t_us(1)) / 1e6;

    % Sequence number (bytes 3-6)
    data.seq = uint32(frames(:,3)) * 2^24 + uint32(frames(:,4)) * 2^16 + ...
               uint32(frames(:,5)) * 2^8  + uint32(frames(:,6));

    % XYZ (bytes 11-19, 20-bit left-justified)
    uX = bitshift(uint32(frames(:,11)),16) + bitshift(uint32(frames(:,12)),8) + uint32(frames(:,13));
    uY = bitshift(uint32(frames(:,14)),16) + bitshift(uint32(frames(:,15)),8) + uint32(frames(:,16));
    uZ = bitshift(uint32(frames(:,17)),16) + bitshift(uint32(frames(:,18)),8) + uint32(frames(:,19));

    uX = bitshift(uX, -4);
    uY = bitshift(uY, -4);
    uZ = bitshift(uZ, -4);

    % 20-bit two's complement
    X_raw = double(uX); X_raw(uX >= 2^19) = X_raw(uX >= 2^19) - 2^20;
    Y_raw = double(uY); Y_raw(uY >= 2^19) = Y_raw(uY >= 2^19) - 2^20;
    Z_raw = double(uZ); Z_raw(uZ >= 2^19) = Z_raw(uZ >= 2^19) - 2^20;

    LSB_PER_G = 1 / 51.2e-6;
    data.X = X_raw / LSB_PER_G;
    data.Y = Y_raw / LSB_PER_G;
    data.Z = Z_raw / LSB_PER_G;

    % TTL (flags byte 20, bit 0)
    data.TTL = double(bitand(frames(:,20), 1) ~= 0);

    % Sampling rate
    if length(data.t) > 1
        data.Fs = 1 / median(diff(data.t));
    else
        data.Fs = 0;
    end

    data.nFrames = size(frames, 1);
end

function CRC_TABLE = build_crc_table()
    CRC_TABLE = zeros(256, 1, 'uint16');
    for i = 0:255
        crc = uint16(i);
        for j = 1:8
            if bitand(crc, 1)
                crc = bitxor(bitshift(crc, -1), uint16(hex2dec('A001')));
            else
                crc = bitshift(crc, -1);
            end
        end
        CRC_TABLE(i+1) = crc;
    end
end