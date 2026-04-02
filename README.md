# Reproducibility package for "A Statistical Framework for Signal Extraction in Noisy Ionospheric Observations from the International Space Station"

This repository is solely for reproducing the data processing and figure generation for the associated paper from provided raw inputs. This is not a 
general-purpose pipeline and is not supported for other datasets.
Run everything via run_all.R (other scripts assume variables loaded by run_all.R and are not standalone).

## Abstract

The Electric Propulsion Electrostatic Analyzer Experiment (ÈPÈE) is a compact ion energy bandpass filter deployed on the International Space Station (ISS) in 
March 2023 and providing continuous measurements through April 2024. This period coincides with the Solar Cycle 25 maximum, capturing unique observations of 
solar activity extremes in the mid- to low-latitude regions of the topside ionosphere.  From these in situ spectra we derive plasma parameters that inform 
space-weather impacts on satellite navigation and radio communication. We present a statistical processing pipeline for ÈPÈE that (i) estimates the instrument 
noise floor, (ii) accounts for irregular temporal sampling, and (iii) extracts ionospheric signals. Rather than discarding noisy data, the method learns a 
baseline noise model and fits the measurement surface using a scaled Vecchia Gaussian process approximation, recovering values typically rejected by 
thresholding. The resulting products increase data coverage and enable noise-assisted monitoring of ionospheric variability.

## Quick start

First, unzip the .zip files in data/epee and data/LLA.

### From the project root (folder with run_all.R)

```
install.packages("renv")
renv::restore()        # installs the exact package versions
source("run_all.R")    # builds data products, figures, and tables
```

### Command line alternative

```
R -q -e "install.packages('renv'); renv::restore(); source('run_all.R')"
```

## Repository structure

```
.
├─ R/
│  ├─ process_funcs.R            # helpers for reading, wrangling, fitting
│  ├─ process.R                  # parallel day-by-day processing (writes logs/plots/results)
│  ├─ post_processing_funcs.R    # helpers for figures/tables
│  └─ post_processing.R          # generates paper figures/tables
├─ data/
│  ├─ epee/                      # raw EPEE inputs
│  ├─ FPMU/                      # raw FPMU inputs
│  ├─ LLA/                       # raw lat/lon/alt inputs
├─ figures/
│  ├─ paper/                     # final, publication-ready figures
│  └─ process/                   # diagnostic plots from processing
├─ logfiles/                     # per-day run logs (created by process.R)
├─ results/                      # smoothed, noise-removed outputs (e.g., *_df.RData)
├─ renv/                         # project-local package library (managed by renv)
├─ renv.lock                     # pinned package versions for reproducibility
└─ run_all.R                     # entrypoint script (loads pkgs, sets paths, sources R/ code)
```

### Expected inputs

- Raw files are under data/epee, data/FPMU, data/LLA, and data/eclipse as shown above.

- Paths are relative to the project root (handled with {here}); no absolute paths are required.


### Outputs

Running run_all.R creates (folders are auto-created if missing):

- Figures
    - figures/paper/ — final figures used in the paper
    - figures/process/ — diagnostic plots created during fitting
- Results 
    - results/ — RData of smoothed, noise-removed current data
- Logs
    - logfiles/ — per-day logs with progress messages


## How to run (important)

Run only via run_all.R.
The other scripts (R/process.R, R/post_processing.R) assume packages, paths, and helper functions have already been loaded by run_all.R. 
Executing them directly may error.


## Software environment

R version 4.5.1 (R Core Team 2025). Packages pinned via renv.lock. If you change packages, run renv::snapshot() and commit. 

Note (Mac only): The parallel code uses a PSOCK cluster (makeCluster()) and has been tested on macOS. Windows/Linux users may need 
to adjust the cluster setup (e.g., FORK on Linux/macOS or additional clusterExport() calls on Windows). This repository is intended to 
reproduce results on macOS systems.

## Copyright and License

Triad National Security, LLC owns the copyright for this code and related data. Please see LICENSE for copyright details.
