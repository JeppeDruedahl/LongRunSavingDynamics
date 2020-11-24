clear all
set more off

//////////////
// 1. setup //
//////////////

set obs 14
gen time = _n + 95
label var time "Years since inheritance"
label define time 0 "96"
forval t = 96/109 {
	label define time `t' "`=`t'-100'", add
}
label values time time



///////////////////////////////////
// 2. data for simulated results //
///////////////////////////////////
	// 2.a Full sample //
	tsset time
	gen zerovar = 0
	foreach model in $old_models $new_models {
		foreach suffix in "" "_unexp" {
			if strpos("$new_models ", "`model' ") {
				local varstobuild A B NW
			}
			else {
				local varstobuild A
			}
			foreach age in all {
				foreach var in `varstobuild' {
					quietly {
						gen `var'_`model'`suffix'_`age' = .
						est use data/local/estimates/simulated/NW_`model'`suffix'_`age'
						scalar toscale = _b[100.timefid]
						est use data/local/estimates/simulated/`var'_`model'`suffix'_`age'
						
						replace `var'_`model'`suffix'_`age' = 0 if time==99
						foreach t of numlist 96/98 {
							replace `var'_`model'`suffix'_`age' = _b[`t'.timefid]/toscale if time==`t'
						}
						foreach t of numlist 100/109 {
							replace `var'_`model'`suffix'_`age' = _b[`t'.timefid]/toscale if (time-1==`t')
						}
						replace `var'_`model'`suffix'_`age' = . if time==100
					}
				}
			
			}
		}
	}


////////////////////////////////////
// 3. data from empirical results //
////////////////////////////////////

 // 3.a Full sample //
	foreach var in networth liqworth {	
		gen beta_`var' = 0
		gen se_`var' 	 = 0
		gen cilo_`var' = 0
		gen ciup_`var' = 0
		est use data/local/estimates/empirical/networth_2_main
		scalar toscale = _b[101.inhgroup_all]
		est use data/local/estimates/empirical/`var'_2_main 
		foreach t of numlist 96/98 100/109 {
			quietly {
				replace beta_`var' = _b[`t'.inhgroup_all]/toscale if time==`t'
				replace se_`var'   = _se[`t'.inhgroup_all]/toscale if time==`t'
				replace ciup_`var' = (_b[`t'.inhgroup_all] + 1.96*_se[`t'.inhgroup_all])/toscale if time==`t'
				replace cilo_`var' = (_b[`t'.inhgroup_all] - 1.96*_se[`t'.inhgroup_all])/toscale if time==`t'
			}
		}
	}
*

///////////////
// 4. graphs //
///////////////

/* Buffer stock, rhos*/
twoway 	(connected A_bs_ez_LCP_rho_1_all time, color(red) ms(+)) ///
		(connected A_bs_ez_LCP_rho_2_all time, color(red) ms(X)) ///
		(line A_bs_ez_LCP_rho_3_all time, color(red) lp(-)) ///
		(line A_bs_ez_LCP_rho_4_all time, color(red) lp(longdash_dot) ) ///
		(scatter beta_networth time, mc(black)) ///
		(rcap ciup_networth cilo_networth time, lc(black)) ///
		, graphr(color(white) margin(tiny)) xlabel(96/109, value) ylabel(-0.2(0.2)1.2) yscale(range(-0.3 1.3)) ///
	xline(100, lc(black) lp(-)) yline(0, lw(medthick) lc(black)) ytitle("Wealth change ratio") ///
	title("Long-run saving dynamics (LRD)", color(black)) legend(r(2) symplacement(left) ///
		region(color(white)) ///
		order(	5 "Empirical estimates" - - - "" - "" ///
				- "Targeting only LCP:" 1 "Risk aversion of1.5" 2 "Risk aversion of2" 3 "Risk aversion of4" 4 "Risk aversion of6" ///
				))
	graph save data/local/graphs/gph/bs_irf_rhos.gph, replace

/* Buffer stock, matched*/
twoway 	(connected A_bs_ez_LCP_IRF_all time, color(blue) lp(-) ms(Oh)) ///
		(connected A_bs_ez_LCP_IRF_per_all time, color(blue) lp(-) ms(Th)) ///
		(connected A_bs_ez_LCP_IRF_var_all time, color(blue) lp(-) ms(Sh)) ///
		(scatter beta_networth time, mc(black)) ///
		(rcap ciup_networth cilo_networth time, lc(black)) ///
		, graphr(color(white) margin(tiny)) xlabel(96/109, value) ylabel(-0.2(0.2)1.2) yscale(range(-0.3 1.3)) ///
	xline(100, lc(black) lp(-)) yline(0, lw(medthick) lc(black)) ytitle("Wealth change ratio") ///
	title("Long-run saving dynamics (LRD)", color(black)) legend(r(2) symplacement(left) ///
		region(color(white)) ///
		order( 4 "Empirical estimates" - - "" - "" ///
				- "Targeting LCP and LRD:" 1 "Free risk aversion" 2 "Free beliefs" 3 "Free income risk" ///
				))
	graph save data/local/graphs/gph/bs_irf_matched.gph, replace

/* kaplan violante*/
twoway 	(connected NW_kv_ez_LCP_rho_2_all time, color(red) ms(X)) ///
		(connected NW_kv_ez_LCP_IRF_all time, color(blue) lp(-) ms(Oh)) ///
		(scatter beta_networth time, mc(black)) ///
		(rcap ciup_networth cilo_networth time, lc(black)) ///
		, graphr(color(white) margin(tiny)) xlabel(96/109, value) ylabel(-0.2(0.2)1.2) yscale(range(-0.3 1.3)) ///
	xline(100, lc(black) lp(-)) yline(0, lw(medthick) lc(black)) ytitle("Wealth change ratio") ///
	title("Long-run saving dynamics (LRD), net worth", color(black)) legend(r(1) symplacement(left) ///
		region(color(white)) ///
		order(1 "Targeting only LCP ({&rho}=2)" 2 "Targeting LCP and LRD (net worth)" 3 "Empirical estimates"))
	graph save data/local/graphs/gph/kv_irf_NW.gph, replace
	
twoway 	(connected A_kv_ez_LCP_rho_2_all time, color(red) ms(X)) ///
		(connected A_kv_ez_LCP_IRF_all time, color(blue) lp(-) ms(Oh)) ///
		(scatter beta_liqworth time, mc(black)) ///
		(rcap ciup_liqworth cilo_liqworth time, lc(black)) ///
		, graphr(color(white) margin(tiny)) xlabel(96/109, value) ylabel(-0.2(0.2)1.2) yscale(range(-0.3 1.3)) ///
	xline(100, lc(black) lp(-)) yline(0, lw(medthick) lc(black)) ytitle("Wealth change ratio") ///
	title("Long-run saving dynamics (LRD), liquid worth", color(black)) legend(r(1) symplacement(left) ///
		region(color(white)) ///
		order(1 "Targeting only LCP ({&rho}=2)" 2 "Targeting LCP and LRD (net worth)" 3 "Empirical estimates"))
	graph save data/local/graphs/gph/kv_irf_A.gph, replace

grc1leg data/local/graphs/gph/kv_irf_NW.gph data/local/graphs/gph/kv_irf_A.gph, ///
	graphr(color(white))
graph display, xsize(10) ysize(4)
graph export data/local/graphs/kv_fits.pdf, replace

/* Appendix - liquid VS wealth shock */
twoway 	(connected A_bs_ez_LCP_rho_2_all time, color(red) ms(X)) ///
		(line A_bs_ez_LCP_rho_2_unexp_all time, lp(-) color(red)) ///
		(scatter beta_networth time, mc(black)) ///
		(rcap ciup_networth cilo_networth time, lc(black)) ///
		, graphr(color(white) margin(tiny)) xlabel(96/109, value) ylabel(-0.2(0.2)1.2) yscale(range(-0.3 1.3)) ///
	xline(100, lc(black) lp(-)) yline(0, lw(medthick) lc(black)) ytitle("Wealth change ratio") ///
	legend(r(1) symplacement(left) ///
		region(color(white)) ///
		order(1 "Exp. inheritance" 2 "Unexp. inheritance" 3 "Empirical estimates"))
graph export data/local/graphs/appendix_bs_expected.pdf, replace

/* Appendix - liquid VS wealth shock */
twoway 	(connected NW_kv_ez_LCP_rho_2_all time, color(red) ms(X)) ///
		(line NW_kv_ez_LCP_rho_2_unexp_all time, lp(-) color(red)) ///
		(scatter beta_networth time, mc(black)) ///
		(rcap ciup_networth cilo_networth time, lc(black)) ///
		, graphr(color(white) margin(tiny)) xlabel(96/109, value) ylabel(-0.2(0.2)1.2) yscale(range(-0.3 1.3)) ///
	xline(100, lc(black) lp(-)) yline(0, lw(medthick) lc(black)) ytitle("Wealth change ratio") ///
	legend(r(1) symplacement(left) ///
		region(color(white)) ///
		order(1 "Exp. inheritance" 2 "Unexp. inheritance" 3 "Empirical estimates"))
graph export data/local/graphs/appendix_kv_expected.pdf, replace
