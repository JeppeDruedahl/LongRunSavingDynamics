# LongRunSavingDynamics

Code for producing all results for "Long-Run Saving Dynamics: Evidence from Unexpected Inheritances", [Druedahl](http://web.econ.ku.dk/druedahl) and [Martinello](http://www.alemartinello.com), *The Review of Economics and Statistics*.

The data used in the paper are drawn from Danish administrative registers and are confidential. However, our access is not unique and others can gain similar access by following a procedure described by Statistics Denmark at [http://www.dst.dk/en/TilSalg/Forskningsservice.aspx](http://www.dst.dk/en/TilSalg/Forskningsservice.aspx). Additional datails are provided in Appendix F in the Supplemental Materials.

The results are reproduced in the following steps:

1. `data/server/master.do` reproduces all stata result files, data, tables and figures that are produced on Danish administrative data. These include the empirical results and the necessary outputs for running the rest of the code (see `data/server/ReadMe.md` for additional details). These files use sensitive administrative data stored in a secure server and are therefore not directly reproducible from this repository alone. These files reproduce all the empirical results in the paper. Crucially, for the model part of the paper, five files need to be exported from the server to a local machine (`income_stds.csv`, `ageatdeath.csv`, `all.csv`, `empirical_networth.csv` and `empirical_liqworth.csv`). These files are already copied in `model/data/*` for convenience.
2. Copy the content of `data/server/toexport/*`, (`income_stds.csv`, `ageatdeath.csv`, `all.csv`, `empirical_networth.csv` and `empirical_liqworth.csv`) to `model/data/*`. This step mimics the process of exporting the necessary moments for estimation from a secure server to a local machine. From this step forward no sensitive data is used, and therefore the code is fully reproducible.
3. Run `model/run.m` to estimate and simulate the structural consumption models and produce some model figures (see `model/ReadMe.md` for additional details).
4. Adjust the path in `data/local/master.do` to point to `code/csv/`. Run `data/local/master.do` to produce the final figures.

All source-files for tables and figures have now been produced. See the list below.

&nbsp;

## Source-files for tables and figures

&nbsp;

### Tables

&nbsp;

**Table 1:**

- `data/local/server/tabs/descriptives.txt`

**Table 2:**

- `data/local/server/tabs/empirics_2.txt`
- `data/local/server/tabs/empirics_0.txt`

**Table 3:**

- `data/local/server/tabs/housing_2.txt`

**Table 4:**

- `data/local/server/tabs/other_2_inc.txt`
- `data/local/server/tabs/other_2_pens.txt`
- `data/local/server/tabs/other_2_hhd.txt`

**Table 5:**

- Does not rely on data (except for income shock standard deviations from `income_stds.csv`)

**Table 6:**

- `model/figs_tabs/table_lines/bs_ez_LCP_rho_2_alphas.txt`
- `model/figs_tabs/table_lines/bs_ez_LCP_rho_3_alphas.txt`
- `model/figs_tabs/table_lines/bs_ez_LCP_rho_4_alphas.txt`
- `model/figs_tabs/table_lines/bs_ez_LCP_IRF_alphas.txt`
- `model/figs_tabs/table_lines/bs_ez_LCP_IRF_var_alphas.txt`

**Table 7:**

- `model/figs_tabs/table_lines/bs_ez_LCP_rho_1_iq_irf_liq.txt`
- `model/figs_tabs/table_lines/bs_ez_LCP_rho_2_iq_irf_liq.txt`
- `model/figs_tabs/table_lines/bs_ez_LCP_rho_3_iq_irf_liq.txt`
- `model/figs_tabs/table_lines/bs_ez_LCP_rho_4_iq_irf_liq.txt`
- `model/figs_tabs/table_lines/bs_ez_LCP_IRF_iq_irf_liq.txt`
- `model/figs_tabs/table_lines/bs_het_LCP_rho_1_iq_irf_liq.txt`
- `model/figs_tabs/table_lines/bs_het_LCP_rho_2_iq_irf_liq.txt`
- `model/figs_tabs/table_lines/bs_het_LCP_rho_3_iq_irf_liq.txt`
- `model/figs_tabs/table_lines/bs_het_LCP_rho_4_iq_irf_liq.txt`
- `model/figs_tabs/table_lines/bs_het_LCP_IRF_iq_irf_liq.txt`
- `model/figs_tabs/table_lines/kv_ez_LCP_rho_1_iq_irf_liq.txt`
- `model/figs_tabs/table_lines/kv_ez_LCP_rho_2_iq_irf_liq.txt`
- `model/figs_tabs/table_lines/kv_ez_LCP_rho_3_iq_irf_liq.txt`
- `model/figs_tabs/table_lines/kv_ez_LCP_rho_4_iq_irf_liq.txt`
- `model/figs_tabs/table_lines/kv_ez_LCP_IRF_iq_irf_liq.txt`

&nbsp;

### Tables - Appendices

&nbsp;

**Table A.1:**

- `data/local/server/tabs/FN_VS_event.txt`

**Table C.1:**

- `data/local/server/tabs/FN_VS_event.txt`

**Table D.1:**

- `model/figs_tabs/table_lines/bs_ez_LCP_sigma_rho_1_noalphas.txt`
- `model/figs_tabs/table_lines/bs_ez_LCP_sigma_rho_2_noalphas.txt`
- `model/figs_tabs/table_lines/bs_ez_LCP_sigma_rho_3_noalphas.txt`
- `model/figs_tabs/table_lines/bs_ez_LCP_sigma_rho_4_noalphas.txt`
- `model/figs_tabs/table_lines/bs_ez_LCP_sigma_rho_5_noalphas.txt`
- `model/figs_tabs/table_lines/bs_ez_LCP_sigma_rho_6_noalphas.txt`
- `model/figs_tabs/table_lines/bs_ez_LCP_sigma_rho_7_noalphas.txt`
- `model/figs_tabs/table_lines/bs_ez_LCP_sigma_rho_8_noalphas.txt`
- `model/figs_tabs/table_lines/bs_ez_LCP_sigma_rho_9_noalphas.txt`
- `model/figs_tabs/table_lines/bs_ez_LCP_sigma_rho_10_noalphas.txt`
- `model/figs_tabs/table_lines/bs_ez_LCP_sigma_rho_11_noalphas.txt`
- `model/figs_tabs/table_lines/bs_ez_LCP_sigma_rho_12_noalphas.txt`
- `model/figs_tabs/table_lines/bs_ez_LCP_sigma_rho_13_noalphas.txt`
- `model/figs_tabs/table_lines/bs_ez_LCP_sigma_rho_14_noalphas.txt`
- `model/figs_tabs/table_lines/bs_ez_LCP_sigma_rho_15_noalphas.txt`
- `model/figs_tabs/table_lines/bs_ez_LCP_sigma_rho_16_noalphas.txt`
- `model/figs_tabs/table_lines/bs_ez_LCP_sigma_rho_17_noalphas.txt`
- `model/figs_tabs/table_lines/bs_ez_LCP_sigma_rho_18_noalphas.txt`

**Table D.2:**

- `model/figs_tabs/table_lines/bs_CRRA_LCP_rho_2_alphas.txt`
- `model/figs_tabs/table_lines/bs_CRRA_LCP_rho_3_alphas.txt`
- `model/figs_tabs/table_lines/bs_CRRA_LCP_rho_4_alphas.txt`
- `model/figs_tabs/table_lines/bs_CRRA_LCP_IRF_var_alphas.txt`

**Table D.3:**

- `model/figs_tabs/table_lines/bs_ez_LCP_IRF_per_rho_1_alphas.txt`
- `model/figs_tabs/table_lines/bs_ez_LCP_IRF_per_rho_2_alphas.txt`
- `model/figs_tabs/table_lines/bs_ez_LCP_IRF_per_rho_3_alphas.txt`
- `model/figs_tabs/table_lines/bs_ez_LCP_IRF_per_alphas.txt`
- `model/figs_tabs/table_lines/bs_ez_LCP_IRF_var_rho_1_alphas.txt`
- `model/figs_tabs/table_lines/bs_ez_LCP_IRF_var_rho_2_alphas.txt`
- `model/figs_tabs/table_lines/bs_ez_LCP_IRF_var_rho_3_alphas.txt`
- `model/figs_tabs/table_lines/bs_ez_LCP_IRF_var_alphas.txt`

**Table D.4:**

- `model/figs_tabs/table_lines/robustness.txt`

**Table F1:**

- Does not rely on data

**Table F2:**

- Does not rely on data

**Table F3:**

- Does not rely on data

**Table G1:**

- `data/local/server/tabs/role_lqconstr_2_rev.txt`

**Table G2:**

- `data/local/server/tabs/app_main.txt`

**Table G3:**

- `data/local/server/tabs/app_main_abs.txt`

**Table G4:**

- `data/local/server/tabs/app_placebo.txt`

**Table G5:**

- `data/local/server/tabs/app_placebo_abs.txt`

**Table G6:**

- `data/local/server/tabs/housing.txt`

**Table G7:**

- `data/local/server/tabs/income_pension.txt`

**Table G8:**

- `data/local/server/tabs/app_household.txt`

&nbsp;

### Figures

&nbsp;

**Figure 1:**

- `data/server/graphs/example/avgs_fe_1.pdf`
- `data/server/graphs/example/diffs_fe_1.pdf`

**Figure 2:**

- `data/server/graphs/networth_2_main.pdf`
- `data/server/graphs/networth_0_main.pdf`
- `data/server/graphs/liqworth_2_main.pdf`
- `data/server/graphs/wealthcomp_2_main.pdf`

**Figure 3:**

- `data/local/graphs/bs_fits_rhos.pdf`
- `data/local/graphs/bs_fits_matched.pdf`

**Figure 4:**

- `data/local/graphs/bs_prec.pdf`

**Figure 5:**

- `data/local/graphs/kv_fits.pdf`

&nbsp;

### Figures - Appendices

&nbsp;

**Figure C.1:**

- `data/server/graphs/example/avgs_nofe_1.pdf`
- `data/server/graphs/example/avgs_fe_1.pdf`
- `data/server/graphs/example/diffs_nofe_1.pdf`
- `data/server/graphs/example/diffs_fe_1.pdf`

**Figure C.2:**

- `data/server/graphs/example/diffs_fe_1.pdf`
- `data/server/graphs/example/diffs_fe_2.pdf`
- `data/server/graphs/example/eventrestr_fe_1.pdf`
- `data/server/graphs/example/eventrestr_fe_2.pdf`
- `data/server/graphs/example/event_fe_1.pdf`
- `data/server/graphs/example/event_fe_2.pdf`

**Figure C.3:**

- `data/server/graphs/example/networth_FN_yearcohort_balance.pdf`
- `data/server/graphs/example/networth_event_yearcohort.pdf`
- `data/server/graphs/example/liqassets_FN_yearcohort_balance.pdf`
- `data/server/graphs/example/liqassets_event_yearcohort.pdf`

**Figure D.1:**

- `model/figs_tabs/bs_ez_LCP_IRF/inh_cycle.png`
- `model/figs_tabs/bs_ez_LCP_IRF/lcp_inh_p_mean.png`

**Figure D.2:**

- `model/figs_tabs/bs_ez_LCP_IRF/fit_P_mean.png`
- `model/figs_tabs/bs_ez_LCP_IRF/fit_H_pos_mean.png`
- `model/figs_tabs/bs_ez_LCP_IRF/fit_h_cond_mean.png`

**Figure D.3:**

- `data/local/appendix_bs_expected.pdf`
- `data/local/appendix_kv_expected.pdf`

**Figure E.1:**

- `model/figs_tabs/test_suite/test_bs_PIH_lcp_compare_C_mean.png`
- `model/figs_tabs/test_suite/test_bs_PIH_lcp_compare_AB_mean.png`

**Figure E.2:**

- `model/figs_tabs/test_suite/test_bs_vfi_negm_lcp_compare_C_mean.png`
- `model/figs_tabs/test_suite/test_bs_vfi_negm_lcp_compare_AB_mean.png`

**Figure E.3:**

- `model/figs_tabs/test_suite/test_bs_ez_compare_AB_t35_ez.png`
- `model/figs_tabs/test_suite/test_bs_zeta_compare_AB_t35_zetas.png`

**Figure E.4:**

- `model/figs_tabs/test_suite/test_bs_grid_compare_AB_t20_grid.png`
- `model/figs_tabs/test_suite/test_bs_grid_compare_AB_t35_grid.png`

**Figure E.5:**

- `model/figs_tabs/test_suite/test_lambda_compare_AB_t20_lambas.png`
- `model/figs_tabs/test_suite/test_lambda_compare_AB_t35_lambas.png`

**Figure E.6:**

- `model/figs_tabs/test_suite/test_grid_compare_AB_t20_grid.png`
- `model/figs_tabs/test_suite/test_grid_compare_AB_t35_grid.png`

**Figure E.7:**

- `model/figs_tabs/test_suite/test_vfi_negm_lcp_compare_C_mean.png`
- `model/figs_tabs/test_suite/test_vfi_negm_lcp_compare_AB_mean.png`

**Figure E.8:**

- `model/figs_tabs/test_suite/test_singles_beta_compare_AB_t35_beta.png`
- `model/figs_tabs/test_suite/test_singles_rho_compare_AB_t35_rho.png`
- `model/figs_tabs/test_suite/test_singles_sigma_compare_AB_t35_sigma.png`
- `model/figs_tabs/test_suite/test_singles_zeta_compare_AB_t35_zeta.png`
- `model/figs_tabs/test_suite/test_singles_kappa_compare_AB_t35_kappa.png`
- `model/figs_tabs/test_suite/test_singles_sigma_psi_compare_AB_t35_sigma_psi.png`

**Figure G.1:**

- `data/server/graphs/robustness_byinher.pdf`
- `data/server/graphs/robustness_byinher_norm.pdf`

**Figure G.2:**

- `data/server/graphs/inheritance_distribution.pdf`

**Figure G.3:**

- `data/server/graphs/distribution_perminc.pdf`
- `data/server/graphs/distribution_networth.pdf`

**Figure G.4:**

- `data/server/graphs/robustness_norm_VS_abs.pdf`

**Figure G.5:**

- `data/server/graphs/robustness_base_VS_inflated_housing.pdf`

**Figure G.6:**

- `data/server/graphs/robustness_networth_individual_VS_hhd.pdf`
- `data/server/graphs/robustness_liqworth_individual_VS_hhd.pdf`

**Figure G.7:**

- `data/server/graphs/robustness_networth_liquid_vs_illiquid.pdf`
