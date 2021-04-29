## Overview

This repository contains code to produce the results from "Fission-fusion dynamics and social associations are structured by daily ranging and den usage in spotted hyenas"

The README describes the scripts used in the analysis and how to implement them. 

## Scripts

#### Main analyses
1. ff_driver.R - Main driver script containing the function that runs all analyses and calls of this function. This is also where you need to specify necessary filepaths. 
2. ff_functions_library.R - Library of functions called by main driver function for running specific analyses and producing plots. This script is sourced by the driver function. 
3. ff_summary_figures.R - This script produces prints summary statistics and outputs plots that appear in the manuscript. 

#### Preprocessing
1. hyena_preprocess_0_extract_GPS_csv_to_R.R - preprocessing step for extracting and formatting GPS data from collar output
2. hyena_preprocess_1_filter_gps.R - preprocessing step for filtering GPS data
3. hyena_preprocess_2_link_gps_and_vedba.R - preprocessing step for linking GPS and accelerometer data

## How to run the analysis
Here are instructions for how to run this analysis. 

To run the analysis you will need to specify four file paths and pass them to the driver function in ff_driver.R:
1. code.directory - filepath to code contained in this repo
2. raw.data.directory - filepath 
3. processed.data.directory - filepath to results folder where processed data will be output.
4. results.directory - filepath to results folder where analyzed data and plots will be output. This folder must contain specific subdirectories to store output for runs with different parameters.

After defining these directories, you are ready to run the analysis. Here are the arguments for the runall() driver function:

runall <- function(R.fusion = 100, 
                   R.fission = 200,
                   n.rands = 100,
                   raw.data.directory, processed.data.directory, results.directory, code.directory,
                   ensure.no.day.matches = T, 
                   verbose = T,
                   preprocess = T, ## The rest are options for skipping steps to reduce runtime
                   extract.ff.events = T,
                   get.ff.features = T,
                   execute.day.randomization = T,
                   get.sync.measures = T){