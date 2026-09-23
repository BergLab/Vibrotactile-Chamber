# Vibrotactile Chamber Analysis

MATLAB pipeline for the mechanical characterisation of a five-position
vibrotactile stimulation platform. Raw tri-axial accelerometer recordings
are parsed, segmented into stimulation bursts and trains, and analysed for
displacement, mechanical delay, rise time, decay time, train dynamics and
spectral purity.

## Dependencies

Developed and tested in MATLAB R2025b. Requires the following MathWorks
toolboxes:

- Signal Processing Toolbox (`butter`, `filtfilt`, `hilbert`, `findpeaks`, `spectrogram`)
- Wavelet Toolbox (`cwt`, used in the spectral-purity analysis)
- Statistics and Machine Learning Toolbox (`iqr`)

## Data

The pipeline expects one binary file per platform position, all located in a
single folder. Files are matched by name in the order `Pos1_Accelerometers.bin`
… `Pos5_Accelerometers.bin`, with a fallback to the physical position labels
`OST_`, `SYD_`, `VEST_`, `NORD_`, `MID_Accelerometers.bin`.

Each position file contains a fixed protocol recorded as one continuous stream:
a calibration block (7 potentiometer steps × 3 carrier frequencies × 5
repetitions = 105 bursts) followed by a train block (3 potentiometer steps ×
4 inter-burst intervals × 5 repetitions × 6 bursts per train = 360 bursts).

The data folder is **not** included in this repository and is selected at
runtime through a dialog.

### Binary file format

Each file is a stream of 23-byte frames written by the ADXL357 acquisition
firmware (±10 g range, ~4 kHz sample rate; the sampling frequency is inferred
from the median timestamp interval). Each frame is laid out as:

```
[A5 5A] [seq:4] [timestamp_us:4] [XYZ:9] [flags:1] [CRC16:2] [55]
```

- Header bytes `A5 5A`, tail byte `55`.
- CRC16 (Modbus, polynomial 0xA001) computed over bytes 3–20; frames that fail
  the check are discarded.
- X, Y, Z are 20-bit two's-complement values, left-justified across three bytes
  each, scaled to g at 51.2 µg/LSB.
- The stimulus TTL is bit 0 of the flags byte (1 = stimulation on).
- The Z axis is the primary vibration axis.

`read_adxl_bin.m` parses a single file into a struct (`t`, `X`, `Y`, `Z`,
`TTL`, `Fs`, `seq`, `nFrames`). `preprocess_chamber_data.m` then returns a
1×5 struct array `VibData`, one entry per position, with fields `.Name`,
`.Raw`, `.Markers`, `.Calibration` and `.Trains`.

## Usage

Open **Main_VibrationChamberAnalysis.m** and run it. A dialog prompts for the
data folder; the script then preprocesses the recordings and runs the full
analysis, bundling all results in a single struct `R`.

```matlab
Folder  = uigetdir;                     % select the data folder
VibData = preprocess_chamber_data(Folder);
R       = run_analysis(VibData);        % all six analyses
```

Passing a second argument enables the per-analysis diagnostic plots:

```matlab
R = run_analysis(VibData, struct('plot_all', true));
```

The script's first section (`SETUP PATH`) adds the `HelperFunctions/` folder to
the MATLAB path. Run the whole file with the Run button (F5), or — if running
section-by-section — run the `SETUP PATH` section once at the start of the
session before the others. Alternatively, add the repository folder to the path
manually.

## Pipeline

**preprocess_chamber_data.m** loads each position file via **read_adxl_bin.m**,
detects every stimulus onset with **detect_all_rising_edges.m**, and segments
the recording into calibration bursts and trains. The parsed and validated
data for all five positions are returned in the struct array `VibData`.

**run_analysis.m** is the wrapper that runs the six analyses on `VibData` and
collects their outputs into the result struct `R`:

- **compute_displacement.m** — double-integrates Z-axis acceleration to
  steady-state peak-to-peak and RMS displacement (µm) per burst.
- **compute_mechanical_delay.m** — latency between the commanded stimulus
  (TTL) and the mechanical onset in the accelerometer signal.
- **compute_rise_time.m** — envelope-based rise time from mechanical onset to
  90 % of the steady-state amplitude.
- **compute_decay_time.m** — decay time constant and settling time after
  stimulus offset.
- **compute_train_decay.m** — train dynamics: per-burst rise, residual
  amplitude across the inter-burst interval, and final settling.
- **compute_spectral_purity.m** — steady-state total harmonic distortion
  (THD) and fundamental power, with CWT/STFT time–frequency views.

Results are accessed through the fields of `R` (`R.displacement`, `R.delay`,
`R.rise`, `R.decay`, `R.train`, `R.spectral`).

## Repository structure

```
Main_VibrationChamberAnalysis.m   entry script
run_analysis.m                    analysis wrapper
HelperFunctions/                  preprocessing and the six analyses
```

`Main_VibrationChamberAnalysis.m` adds `HelperFunctions/` to the MATLAB path
automatically when the script is run (see Usage for section-by-section
execution).

## Authors

R. J. F. Sørensen, M. C. D. Larsen
