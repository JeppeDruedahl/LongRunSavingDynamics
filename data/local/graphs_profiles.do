clear all
set more off
*cd "C:\Dropbox\Projects\inheritance\Buffer_stock\ale+jeppe\code\stata"


*************************
******* DATABUILD *******
*************************
import delimited model/data/all.csv, clear 

keep alder _networth* _liqworth*
rename alder age
rename _* *

merge 1:1 age using data/local/datacustom/profiles
drop _merge

tsset age

/* 
	nr = no retirement
	np = no precautionary 
*/

*local todecompose = ustrregexra("$old_models", "([a-z]|[A-Z]|_)+_nr", "")
*local todecompose = ustrregexra("`todecompose'", "prefered", "preferred")
local todecompose 	bs_ez_LCP_rho_2 bs_ez_LCP_IRF bs_ez_LCP_IRF_per bs_ez_LCP_IRF_var
cap drop prec_*
cap drop lfc_*
foreach model in `todecompose' {
	/*
	gen prec_w_`model' = mean_NW_`model'_nr
	gen lfc_w_`model' = mean_NW_`model' - mean_NW_`model'_nr
	
	gen prec_s_`model' = d.prec_w_`model'
	gen lfc_s_`model' = d.lfc_w_`model'
	*/
	gen prec_w_`model' = (mean_NW_`model'_nr /*- mean_NW_`model'_nvr*/)
	gen lfc_w_`model' = (mean_NW_`model' - mean_NW_`model'_nr)
	gen prec_p_`model' = prec_w_`model'/(mean_NW_`model' - mean_NW_`model'_nvr)
	gen lfc_p_`model' = lfc_w_`model'/(mean_NW_`model' - mean_NW_`model'_nvr)
}

*



**********************************************
/********** Fit of median profiles **********/
**********************************************
* NOTE: Need to run irfs file first */
// Buffer-stock - Rhos
twoway	(connected p50_NW_bs_ez_LCP_rho_1 age, color(red) ms(+)) ///
		(connected p50_NW_bs_ez_LCP_rho_2 age, color(red) ms(X)) ///
		(line p50_NW_bs_ez_LCP_rho_3 age, color(red) lp(-)) ///
		(line p50_NW_bs_ez_LCP_rho_4 age, color(red) lp(longdash_dot) ) ///
		(connected networth_p50 age, color(black)) ///
	, graphr(color(white) margin(tiny)) xtitle("Age") ytitle("Years of permanent income") ///
	xlabel(25(5)60) /// 
	title("Life-cycle wealth profile (LCP)", color(black)) legend(r(2) symplacement(left) ///
		region(color(white)) ///
		order(	- "" 5 "Empirical estimates" - - "" - "" ///
				- "Targeting only LCP:" 1 "Risk aversion of 1.5" 2 "Risk aversion of 2" 3 "Risk aversion of 4" 4 "Risk aversion of 6" ///
				))
graph save data/local/graphs/gph/bs_lcp_rhos, replace

grc1leg data/local/graphs/gph/bs_lcp_rhos.gph data/local/graphs/gph/bs_irf_rhos.gph, ///
	graphr(color(white)) //legendfrom(data/local/graphs/gph/bs_irf.gph)
graph display, xsize(10) ysize(4)
graph export data/local/graphs/bs_fits_rhos.pdf, replace

// Buffer-stock matched
twoway	(connected p50_NW_bs_ez_LCP_IRF age, color(blue) lp(-) ms(Oh)) ///
		(connected p50_NW_bs_ez_LCP_IRF_per age,color(blue) lp(-) ms(Th)) ///
		(connected p50_NW_bs_ez_LCP_IRF_var age, color(blue) lp(-) ms(Sh)) ///
		(connected networth_p50 age, color(black)) ///
	, graphr(color(white) margin(tiny)) xtitle("Age") ytitle("Years of permanent income") ///
	xlabel(25(5)60) /// 
	title("Life-cycle wealth profile (LCP)", color(black)) legend(r(2) symplacement(left) ///
		region(color(white)) ///
		order(	- "" 4 "Empirical estimates" - - "" - "" ///
				- "Targeting both LCP and LRD:" 1 "Free risk aversion" 3 "Free income risk" 2 "Free perceived income risk" ///
				))
graph save data/local/graphs/gph/bs_lcp_matched, replace

grc1leg data/local/graphs/gph/bs_lcp_matched.gph data/local/graphs/gph/bs_irf_matched.gph, ///
	graphr(color(white)) //legendfrom(data/local/graphs/gph/bs_irf.gph)
graph display, xsize(10) ysize(4)
graph export data/local/graphs/bs_fits_matched.pdf, replace


***************************************************
/***** Precautionary/lifecycle decomposition *****/
***************************************************

/* Buffer-stock */
twoway	(connected prec_w_bs_ez_LCP_rho_2 age, color(red) ms(X)) ///
		(connected prec_w_bs_ez_LCP_IRF age, color(blue) lp(-) ms(Oh)) ///
		(connected prec_w_bs_ez_LCP_IRF_per age, color(blue) lp(-) ms(Th)) ///
		(connected prec_w_bs_ez_LCP_IRF_var age, color(blue) lp(-) ms(Sh)) ///
	if inrange(age,30,60), graphr(color(white) margin(tiny)) xtitle("Age") ///
	xlabel(30(5)60) ylabel(-0(.1)0.9) yline(0, lw(medthick) lc(black)) ytitle("Years of permanent income") ///
	title("Precautionary savings", color(black)) legend(r(2) symplacement(left) ///
		region(color(white)) ///
		order(	- "Fitting only LCP:" 1 "Risk aversion of 2" - "" - "" ///
				- "Fitting LCP and LRD:" 2 "Free risk aversion" 3 "Free perceived income risk" 4 "Free income risk"))
	graph save data/local/graphs/gph/bs_prec.gph, replace
	
twoway	(connected lfc_w_bs_ez_LCP_rho_2 age, color(red) ms(X)) ///
		(connected lfc_w_bs_ez_LCP_IRF age, color(blue) lp(-) ms(Oh)) ///
		(connected lfc_w_bs_ez_LCP_IRF_per age, color(blue) lp(-) ms(Th)) ///
		(connected lfc_w_bs_ez_LCP_IRF_var age, color(blue) lp(-) ms(Sh)) ///
	if inrange(age,30,60), graphr(color(white) margin(tiny)) xtitle("Age") ///
	xlabel(30(5)60) ylabel(-0.25(.25)2.25) yline(0, lw(medthick) lc(black)) ytitle("Years of permanent income") ///
	title("Residual (life-cycle) savings", color(black) ) legend(off)
	graph save data/local/graphs/gph/bs_lfc.gph, replace
	
grc1leg data/local/graphs/gph/bs_prec.gph data/local/graphs/gph/bs_lfc.gph, ///
	graphr(color(white))
graph display, xsize(10) ysize(4)
graph export data/local/graphs/bs_prec.pdf, replace

/*
/* Kaplan-Violante */
twoway	(connected prec_w_kv_ez_LCP age, color(red) ms(X)) ///
		(connected prec_w_kv_ez_LCP_IRF age, color(blue) lp(-) ms(Oh)) ///
	if inrange(age,30,60), graphr(color(white) margin(tiny)) xtitle("Age") ///
	xlabel(30(5)60) ylabel(-0(.1)1.2) yline(0, lw(medthick) lc(black)) ytitle("Years of permanent income") ///
	title("Precautionary savings", color(black)) legend(r(1) symplacement(left) ///
		region(color(white)) ///
		order(1 "Fitting only LCP" 2 "Fitting LCP and LRD"))
	graph save data/local/graphs/gph/kv_prec.gph, replace
	
twoway	(connected lfc_w_kv_ez_LCP age, color(red) ms(X)) ///
		(connected lfc_w_kv_ez_LCP_IRF age, color(blue) lp(-) ms(Oh)) ///
	if inrange(age,30,60), graphr(color(white) margin(tiny)) xtitle("Age") ///
	xlabel(30(5)60) ylabel(-0.25(.25)2.25) yline(0, lw(medthick) lc(black)) ytitle("Years of permanent income") ///
	title("Residual (lifecycle) savings", color(black) ) legend(off)
	graph save data/local/graphs/gph/kv_lfc.gph, replace
	
grc1leg data/local/graphs/gph/kv_prec.gph data/local/graphs/gph/kv_lfc.gph, ///
	graphr(color(white))
graph display, xsize(10) ysize(4)
graph export data/local/graphs/kv_prec.pdf, replace

