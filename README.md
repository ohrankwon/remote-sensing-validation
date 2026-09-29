# Code and data for "Algorithm evaluation for ocean colour remote sensing"

This archive reproduces every figure and table in the main paper and the
supporting information. All scripts are R and are meant to be run from the
**top level of this directory** (`code-alg-eval-ocean-submittion/`, the folder
containing this README; they use relative paths, so set the working
directory here). 

## Which script makes which figure/table

Run the numbered scripts in order first (they produce the intermediate results
in `output/RData/`), then run any `plot*.R` script. Each plotting script is
named after the paper output it produces.

| Script | Produces | Reads |
|---|---|---|
| `plot1.R` | **Figure 1** (main) | `analysis_random_mse_clust.RData`, `analysis_random_mse_noclust.RData` |
| `plot2-table1-tableS2.R` | **Figure 2** (main), **Table 1** (main), **Table S2** | `analysis_cross_region.RData` |
| `plotS2-tableS1-tableS3.R` | **Figure S2**, **Table S1**, **Table S3** | `analysis_temporal.RData` |
| `plotS3.R` | **Figure S3** | `analysis_cross_region_grid.RData` |
| `plotS4.R` | **Figure S4** | `analysis_random_mae_clust.RData`, `analysis_random_mae_noclust.RData` |
| `plotS5.R` | **Figure S5** | `analysis_combined.RData` |

Figures are written to `output/figure/`. Tables are printed to the console as
LaTeX (via `xtable`) and were pasted into the manuscript by hand; the numbers
are not otherwise post-processed.

**Figure S1** (`output/figure/map_freq.pdf`, the global map of matched
observations) was produced separately and its plotting script is not part of
this archive. The PDF is included for completeness.

## Pipeline

| Script | Purpose | In | Out |
|---|---|---|---|
| `01-cluster.R` | spatiotemporal clustering (T = 7 days, S = 50 km, single linkage cut at height 1) | `data/chl_data_log_source.csv` | `data/chl_data_log_cluster.csv` |
| `02-analysis-random-mse.R` | within-region validation R1–R1000, MSE; runs both Clust and NoClust | `data/chl_data_log_cluster.csv` | `analysis_random_mse_{clust,noclust}.RData` |
| `03-analysis-random-mae.R` | as above with MAE instead of MSE | `data/chl_data_log_cluster.csv` | `analysis_random_mae_{clust,noclust}.RData` |
| `04-analysis-cross-region.R` | cross-region validation E1–E4 | `data/chl_data_log_cluster.csv` | `analysis_cross_region.RData` |
| `05-analysis-cross-region-grid.R` | all 16 train/test zone combinations | `data/chl_data_log_cluster.csv` | `analysis_cross_region_grid.RData` |
| `06-analysis-temporal.R` | temporal validation T1–T7 | `data/chl_data_log_cluster.csv` | `analysis_temporal.RData` |
| `07-analysis-combined.R` | combined cross-region + temporal, scenarios M1/M2 | `data/chl_data_log_cluster.csv` | `analysis_combined.RData` |

`R/functions.R` holds the two helpers shared across scripts: `predict_1nn()`
and `spatiotemporal_cluster()`.

Step 01 only needs to be re-run if you want to rebuild
`data/chl_data_log_cluster.csv`; that file is already included, so scripts
02–07 can be run directly. 

## Data

| File | Description |
|---|---|
| `data/Valente_2022_insitu_chla.csv` | in-situ chlorophyll compilation of Valente et al. (2022); included as the original source of the measurements, not read by any script |
| `data/chl_data_log_source.csv` | the n = 18,564 matched in-situ / OC-CCI Rrs pairs, log10-transformed, with the in-situ sub-dataset label in `source` — **the input to `01-cluster.R`** |
| `data/chl_data_log_cluster.csv` | as above, plus the cluster id (from `01-cluster.R`) — **the input to all analyses** |

Columns `log_x1`–`log_x6` are log10 Rrs at 412, 443, 490, 510, 560 and 665 nm;
`log_y` is log10 chlorophyll.

The satellite/in-situ matching that produced `data/chl_data_log_source.csv`
(same-day pixels within 1 km of each measurement site, from the OC-CCI v6.0
1 km product) was carried out separately.

## Use of AI assistance

Claude (Anthropic) was used in preparing this archive: reorganising and
documenting the analysis scripts for release, and drafting this README.
The study design, statistical methodology, and interpretation of the results
are the authors' own. All AI-assisted code and documentation were reviewed and
verified by the authors.

## Requirements

R (tested on 4.5.2) with: `ggplot2`, `tidyverse`, `tidyr`, `dplyr`, `scales`,
`xtable`, `ggpattern`, `reshape2`, `matrixStats`, `withr`, `geosphere`.
