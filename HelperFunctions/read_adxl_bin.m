function data = read_adxl_bin(filepath)
% READ_ADXL_BIN  Parse ADXL357 binary file (23-byte frames) into struct.
%
%   data = read_adxl_bin(filepath)
%
%   Frame format: [A5 5A] [seq(4)] [t_us(4)] [XYZ(9)] [flags(1)] [CRC(2)] [55]
%                  1  2     3-6       7-10      11-19     20         21-22    23
%
%   Output struct fields:
%     .t    — time vector (seconds, zeroed to first sample)
%     .X    — X acceleration (g)
%     .Y    — Y acceleration (g)
%     .Z    — Z acceleration (g)
%     .TTL  — logical: stimulus ON flag (bit 0 of flags byte)
%     .Fs   — sampling frequency (Hz, from median dt)
%     .seq  — sequence numbers (uint32)
%     .nFrames — total valid frames

    data = [];

    if ~isfile(filepath)
        warning('File not found: %s', filepath);
        return;
    end

    fid = fopen(filepath, 'rb');
    raw = fread(fid, Inf, '*uint8');
    fclose(fid);

    %% 1. Find candidate frames by header/tail bytes
    FRAME_LEN = 23;
    nBytes = numel(raw);

    if nBytes < FRAME_LEN
        warning('File too small: %s', filepath);
        return;
    end

    % Vectorized search for A5 5A ... 55 pattern
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

    % Build frame matrix (nFrames × 23)
    idx_mat = cand(:) + (0:22);
    frames = reshape(raw(idx_mat.'), FRAME_LEN, []).';
    fprintf('  %s: %d candidate frames', filepath, size(frames,1));

    %% 2. CRC16 validation (Modbus, poly 0xA001)
    CRC_TABLE = build_crc_table();

    n = size(frames, 1);
    crc_calc = repmat(uint16(65535), n, 1);

    for col = 3:20   % bytes 3–20: seq, timestamp, xyz, flags
        byte = uint16(frames(:, col));
        lut_idx = bitxor(bitand(crc_calc, uint16(255)), byte) + 1;
        crc_calc = bitxor(bitshift(crc_calc, -8), CRC_TABLE(lut_idx));
    end

    crc_file = bitor(bitshift(uint16(frames(:,21)), 8), uint16(frames(:,22)));
    valid = (crc_calc == crc_file);
    frames = frames(valid, :);
    fprintf(' → %d valid (%.1f%%)\n', size(frames,1), 100*sum(valid)/n);

    if isempty(frames)
        warning('No CRC-valid frames in: %s', filepath);
        return;
    end

    %% 3. Extract timestamp (bytes 7–10, big-endian uint32 µs)
    t_us = double(frames(:,7)) * 2^24 + double(frames(:,8)) * 2^16 + ...
           double(frames(:,9)) * 2^8  + double(frames(:,10));
    data.t = (t_us - t_us(1)) / 1e6;

    %% 4. Extract sequence number (bytes 3–6, big-endian uint32)
    data.seq = uint32(frames(:,3)) * 2^24 + uint32(frames(:,4)) * 2^16 + ...
               uint32(frames(:,5)) * 2^8  + uint32(frames(:,6));

    %% 5. Extract XYZ (bytes 11–19, 20-bit left-justified in 3 bytes)
    % Combine 3 bytes, right-shift by 4 to get 20-bit value
    uX = bitshift(uint32(frames(:,11)),16) + bitshift(uint32(frames(:,12)),8) + uint32(frames(:,13));
    uY = bitshift(uint32(frames(:,14)),16) + bitshift(uint32(frames(:,15)),8) + uint32(frames(:,16));
    uZ = bitshift(uint32(frames(:,17)),16) + bitshift(uint32(frames(:,18)),8) + uint32(frames(:,19));

    uX = bitshift(uX, -4);
    uY = bitshift(uY, -4);
    uZ = bitshift(uZ, -4);

    % 20-bit two's complement → signed double
    X_raw = double(uX); X_raw(uX >= 2^19) = X_raw(uX >= 2^19) - 2^20;
    Y_raw = double(uY); Y_raw(uY >= 2^19) = Y_raw(uY >= 2^19) - 2^20;
    Z_raw = double(uZ); Z_raw(uZ >= 2^19) = Z_raw(uZ >= 2^19) - 2^20;

    % Scale to g (ADXL357 ±10g: 51.2 µg/LSB)
    LSB_PER_G = 1 / 51.2e-6;   % = 19531.25
    data.X = X_raw / LSB_PER_G;
    data.Y = Y_raw / LSB_PER_G;
    data.Z = Z_raw / LSB_PER_G;

    %% 6. Extract TTL (flags byte 20, bit 0 = FLAG_STIM_ON)
    data.TTL = double(bitand(frames(:,20), 1) ~= 0);

    %% 7. Sampling rate
    if length(data.t) > 1
        data.Fs = 1 / median(diff(data.t));
    else
        data.Fs = 0;
    end

    data.nFrames = size(frames, 1);
end


function CRC_TABLE = build_crc_table()
% Build CRC16 Modbus lookup table (poly 0xA001)
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