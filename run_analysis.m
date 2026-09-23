function R = run_analysis(VibData, opts)
% RUN_ANALYSIS  Run all vibrotactile platform characterization analyses.
%
%   R = run_analysis(VibData)
%   R = run_analysis(VibData, opts)
%
%   Calls all analysis functions, bundles results into a single struct R.
%   Plot generation in each function is controlled via the opts struct.
%
%   Options (all default to false):
%     opts.plot_displacement   — plot displacement scatter/summary
%     opts.plot_delay          — plot mechanical delay diagnostics
%     opts.plot_rise           — plot rise time diagnostics
%     opts.plot_decay          — plot decay time diagnostics
%     opts.plot_train          — plot train decay diagnostics
%     opts.plot_spectral       — plot spectral purity (CWT/STFT/THD)
%     opts.plot_all            — shortcut: set ALL plot flags to true
%
%   Output struct R:
%     R.displacement   — from compute_displacement
%     R.delay          — from compute_mechanical_delay
%     R.rise           — from compute_rise_time
%     R.decay          — from compute_decay_time
%     R.train          — from compute_train_decay
%     R.spectral       — from compute_spectral_purity (position 1 only)
%     R.opts           — the options used
%     R.timestamp      — datestr of when analysis was run
%     R.nPos           — number of positions in VibData
%
%   Example:
%     VibData = preprocess_chamber_data(Folder);
%
%     % Run all analyses, no plots
%     R = run_analysis(VibData);
%
%     % Run all analyses with plots
%     R = run_analysis(VibData, struct('plot_all', true));
%
%     % Run with only decay plots
%     R = run_analysis(VibData, struct('plot_decay', true));

    %% === PARSE OPTIONS ===
    if nargin < 2, opts = struct(); end

    % Default all plot flags to false
    default_opts = struct( ...
        'plot_displacement', false, ...
        'plot_delay',        false, ...
        'plot_rise',         false, ...
        'plot_decay',        false, ...
        'plot_train',        false, ...
        'plot_spectral',     false, ...
        'plot_all',          false  ...
    );

    % Merge user opts over defaults
    fnames = fieldnames(default_opts);
    for i = 1:numel(fnames)
        if ~isfield(opts, fnames{i})
            opts.(fnames{i}) = default_opts.(fnames{i});
        end
    end

    % plot_all overrides individual flags
    if opts.plot_all
        opts.plot_displacement = true;
        opts.plot_delay        = true;
        opts.plot_rise         = true;
        opts.plot_decay        = true;
        opts.plot_train        = true;
        opts.plot_spectral     = true;
    end

    %% === RUN ANALYSES ===
    fprintf('\n');
    fprintf('##############################################################\n');
    fprintf('  VIBROTACTILE PLATFORM CHARACTERIZATION — FULL ANALYSIS\n');
    fprintf('  %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf('##############################################################\n\n');

    % 1. Displacement
    fprintf('>>> [1/6] Displacement calibration...\n');
    R.displacement = compute_displacement(VibData, opts.plot_displacement);

    % 2. Mechanical delay
    fprintf('\n>>> [2/6] Mechanical delay...\n');
    R.delay = compute_mechanical_delay(VibData, opts.plot_delay);

    % 3. Rise time
    fprintf('\n>>> [3/6] Rise time...\n');
    R.rise = compute_rise_time(VibData, opts.plot_rise);

    % 4. Decay time
    fprintf('\n>>> [4/6] Decay time...\n');
    R.decay = compute_decay_time(VibData, opts.plot_decay);

    % 5. Train decay
    fprintf('\n>>> [5/6] Train decay...\n');
    R.train = compute_train_decay(VibData, opts.plot_train);

    % 6. Spectral purity
    fprintf('\n>>> [6/6] Spectral purity...\n');
    R.spectral = compute_spectral_purity(VibData, opts.plot_spectral);

    %% === METADATA ===
    R.opts      = opts;
    R.timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    R.nPos      = numel(VibData);

    fprintf('\n##############################################################\n');
    fprintf('  ANALYSIS COMPLETE — %s\n', R.timestamp);
    fprintf('  Results bundled in struct R with fields:\n');
    fprintf('    .displacement  .delay  .rise  .decay  .train  .spectral\n');
    fprintf('##############################################################\n\n');
end