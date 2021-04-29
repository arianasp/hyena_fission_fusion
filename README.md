## Overview

This repository contains code to produce the results from "Fission-fusion dynamics and social associations are structured by daily ranging and den usage in spotted hyenas"

The README describes the scripts used in the analysis and how to implement them. 

## Scripts

#### Main analyses
1. ff_driver.R - Main driver script containing the `runall()` function that runs all analyses. The `runall()` function is called from this script as well and requires as arguments the input and output directory file paths. 
2. ff_summary_figures.R - This script produces prints summary statistics and outputs plots that appear in the manuscript. It is sourced from the `generate_figures()` function, which takes the directory output of `runall()`
3. ff_functions_library.R - Library of functions called by `runall()` and `generate figures()` that run specific analyses and produce plots. This script is sourced by `runall()`. 

#### Preprocessing
These steps are run in `runall()`
1. hyena_preprocess_0_extract_GPS_csv_to_R.R - preprocessing step for extracting and formatting GPS data from collar output.
2. hyena_preprocess_1_filter_gps.R - preprocessing step for filtering GPS data
3. hyena_preprocess_2_link_gps_and_vedba.R - preprocessing step for linking GPS and accelerometer data

## How to run the analysis
Here are instructions for how to run this analysis. 

To run the analysis you will need to specify four file paths and pass them to `runall()` in ff_driver.R:
1. code.directory - filepath to code contained in this repo
2. raw.data.directory - filepath where raw data are located
3. processed.data.directory - filepath to results folder where processed data will be output.
4. results.directory - filepath to results folder where analyzed data and plots will be output. This folder must contain specific subdirectories to store output for runs with different parameters.

After defining these directories, you are ready to run the analysis. The `runall()` function provides options for altering analysis parameters, specifying directories, and skipping steps to reduce runtime. Preprocessing steps are invariant to the parameter values, so preprocessing only needs to be run once. Many of these steps are computationally expensive, especially preprocessing and the reference model.

The `runall()` function saves the analyzed data and outputs a vector of two directories (`data.outdir` and `plots.outdir`) specific to the parameter values used to run the anlysis. These are passed to the `generate_figures()` function. 

Here are specific details about the arguments for the runall() driver function:

```
runall <- function(
  ### Run parameters
  R.fusion = 100,  # threshold proximity for triggering a fission-fusion event
  R.fission = 200, # threshold proximity defining start and end of ff event
  n.rands = 100, # number of iterations of the reference model to run
  ensure.no.day.matches = T, # should randomizations exclude randomly reproduced original data?
  
  ### input/output directories
  raw.data.directory, # directory of raw data, containing metadata, gps data, acc data, den locations, hyena photograph, satllite map
  processed.data.directory, # directory for output of preprocessing steps
  results.directory, # directory for results, including subdirectories for each parameter combination and data/ and plots/ for data output
  code.directory, # directory where code is located 
  
  ### Run options
  verbose = T, # Should progress be printed to console?
  preprocess = T, # Do you want to run preprocessing and overwrite output? (preprocessing steps are identical regardless of parameter values) 
  extract.ff.events = T, # Do you want to extract fission-fusion events and overwrite output?
  get.ff.features = T, # Do you want to extract features of ff events and overwrite output?
  execute.day.randomization = T, # Do you want to run the reference model and overwrite output? 
  get.sync.measures = T # Do you want to calculate synchrony measures?
){
                   
```

### Full analysis code
Here are the funciton calls that produce the full analysis, copied from the ff_driver.R script. 

```
#-----------------------------------RUN ME-------------------------------------------------
set.seed(43410)
print('--------------------------- MAIN RESULTS ---------------------------------')
output.dirs <- runall(ensure.no.day.matches = T, R.fusion = 100, R.fission = 200, n.rands = 100, 
                      raw.data.directory, processed.data.directory, results.directory, code.directory, preprocess = T,
                      execute.day.randomization = T, extract.ff.events = T, get.sync.measures = T, get.ff.features = T)
generate_figures(output.dirs[1], output.dirs[2], code.directory)

print('--------------------------- CHECK SMALLER THRESHOLD ---------------------------------')
output.dirs <- runall(ensure.no.day.matches = T, R.fusion = 50, R.fission = 100, n.rands = 100, 
                      raw.data.directory, processed.data.directory, results.directory, code.directory, preprocess = F,
                      execute.day.randomization = T, extract.ff.events = T, get.sync.measures = T, get.ff.features = T)
generate_figures(output.dirs[1], output.dirs[2], code.directory)
print('--------------------------- CHECK LARGER THRESHOLD ---------------------------------')
output.dirs <- runall(ensure.no.day.matches = T, R.fusion = 200, R.fission = 300, n.rands = 100, 
                      raw.data.directory, processed.data.directory, results.directory, code.directory, preprocess = F,
                      execute.day.randomization = T, extract.ff.events = T, get.sync.measures = T, get.ff.features = T)
generate_figures(output.dirs[1], output.dirs[2], code.directory)
```