# run_all.R — rebuild everything from raw → figures + results
## 1) pin project root for here::here()
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
here::i_am("run_all.R")

## 2) activate renv using the project root (not the current wd)
act <- file.path(here::here(), "renv", "activate.R")
if (file.exists(act)) source(act)


# Load packages (no installs here; rely on renv::restore())
pkgs <- c("dplyr","ggplot2","stringr","zoo", "factoextra", "patchwork", "viridis", "fuzzyjoin",
          "lubridate", "mgcv","data.table","latex2exp","plotly","here", "dplyr", "magrittr", 
          "purrr", "reshape2", "ncdf4", "GPvecchia", "GpGp", "foreach", "doParallel", "tidyr")
invisible(lapply(pkgs, library, character.only = TRUE))

# ---- paths (all relative to project root) ----
fpmu_pathname     <- here::here("data", "FPMU")
lla_pathname      <- here::here("data", "LLA")
rawdata_pathname  <- here::here("data", "epee")
sv_pathname       <- here::here("results")
figure_pathname   <- here::here("figures", "paper")

# Make sure output dirs exist (git won't track empty folders)
invisible(lapply(
  list(sv_pathname, figure_pathname, here::here("figures","process"), here::here("logfiles")),
  dir.create, recursive = TRUE, showWarnings = FALSE
))


# ---- source functions
source(here::here("R", "process_funcs.R"))
source(here::here("R", "post_processing_funcs.R"))

#source code
source(here::here("R", "process.R"))
source(here::here("R", "post_processing.R"))
