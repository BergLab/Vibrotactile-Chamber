function [onset_indices, onset_times] = detect_all_rising_edges(t, ttl)
% DETECT_ALL_RISING_EDGES  Find every 0→1 transition in a logic signal.
%
%   [idx, times] = detect_all_rising_edges(t, ttl)
%
%   No filtering by gap duration — returns every single rising edge.
%   Use this as the foundation; segmentation logic lives upstream.

    d_ttl = diff(double(ttl));
    onset_indices = find(d_ttl == 1) + 1;   % +1: diff(i) reflects change at i+1
    onset_times   = t(onset_indices);
end