# ReadMe - model

## Requirements

The following software is used:

1. Windows 10 64-bit
2. MATLAB 2018b
3. Microsoft Visual Studio 2017 Community Edition
4. Intel Parallel Studio XE 2018
5. NLopt 2.4.2 (included, else see below)

MATLAB was set to use the Intel C++ compiler (choose this by running `mex -setup C++`).

## Overview

All results are produced by the run_m file. It runs:

1. run_01_main.m: run all calibrations
2. run_02_overview.m: produce figures and tables
3. run_03_test_suite.m: produce validation results for the computational appendix
4. run_04_copy_figs_tabes.m: move the relevant files into the output folder

The underlying .m files are:
 
1. setup.m: functions to set up the model and parameters
2. model.m: functions for solving and simulating the model
3. fun.m: auxilary functions
3. estimate.m: functions for estimating the model.
4. figs.m: functions for producing figures (all produced in figs_tabs/).

The folder structure is:

1. cfuncs/: .cpp files used in MEX
2. csv/: temporarily stored csv-files for simulation
3. data/: data files
4. estimates/: estimates for fast re-production of tables
5. figs_tabs/: figures and tables
5. out/: final output

## NLopt

If you experience trouble with using NLopt consider re-installing as follows

1. go to http://ab-initio.mit.edu/wiki/index.php/NLopt_on_Windows
2. download nlopt-2.4.2-dll64.zip and unzip
3. open "x64 Native Tools Command Prompt for VS 2017" and locate unzipped folder
4. run "lib /def:libnlopt-0.def /machine:x64"
5. copy libnlopt-0.lib to cfuncs\ folder
6. copy libnlopt-0.dll to this folder