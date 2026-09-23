%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  MAIN CHAMBER ANALYSIS — VIBROTACTILE PLATFORM CHARACTERIZATION v3
%  Entry script: select data -> preprocess -> run full analysis.

% Pipeline:
%    preprocess_chamber_data(Folder) -> VibData
%    run_analysis(VibData)           -> R   (all six analyses bundled)
%
%  Authors: R. J. F. Sørensen, M. C. D. Larsen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

%% === 0. SETUP PATH ===
% Add the repository's own subfolders (HelperFunctions) to the MATLAB path.
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir) || contains(scriptDir, tempdir)
    scriptDir = fileparts(which('Main_VibrationChamberAnalysis'));
end
addpath(genpath(scriptDir));

%% === 1. SELECT DATA FOLDER & PREPROCESS === 
% Opens a dialog; select the folder containing the five position files
Folder = uigetdir('', 'Select data folder (five position files)');
if isequal(Folder, 0), error('No folder selected.'); end
VibData = preprocess_chamber_data(Folder);

%% === 2. RUN ALL ANALYSES ===
% Runs displacement, mechanical delay, rise time, decay time, train decay and spectral purity; results are bundled in struct R.
% Pass struct('plot_all', true) as a second argument to also show the per-analysis diagnostic plots.
R = run_analysis(VibData);

%% === 3. ACCESS INDIVIDUAL RESULTS ===
% All results are in R:
%   R.displacement   — displacement per condition
%   R.delay          — mechanical delay per burst
%   R.rise           — rise time per burst
%   R.decay          — decay tau + settle per burst
%   R.train          — train decay analysis
%   R.spectral       — spectral purity / THD

% Example: get mean settle time for 250 Hz, Pot 18
mask = [R.decay.FreqHz] == 250 & [R.decay.PotStep] == 18;
fprintf('250 Hz Pot18 settle: %.1f +/- %.1f ms\n', ...
    mean([R.decay(mask).Settle_Meas_ms], 'omitnan'), ...
    std([R.decay(mask).Settle_Meas_ms], 'omitnan'));