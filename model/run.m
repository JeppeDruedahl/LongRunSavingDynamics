clc;
clear;

%% mex

threads = 32;
model.mex_solve(threads,1);
model.mex_simulate(threads,1);
clc;

%% files

run_01_main
run_02_overview
run_03_test_suite
run_04_copy_figs_tabs