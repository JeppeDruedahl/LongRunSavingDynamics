clear all
set more off
cd "C:\Repositories\LongRunSavingDynamics"

global old_models 	bs_ez_LCP_rho_1 bs_ez_LCP_rho_2 bs_ez_LCP_rho_3 bs_ez_LCP_rho_4 bs_ez_LCP_IRF bs_ez_LCP_IRF_per bs_ez_LCP_IRF_var
					
global new_models 	kv_ez_LCP_rho_2 kv_ez_LCP_IRF

global suffix_list 	_nr _nvr _unexp
									
										
do databuild
do analysis

do graphs_irfs
do graphs_profiles
