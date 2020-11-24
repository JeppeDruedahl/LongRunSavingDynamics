cd "D:\Data\workdata\705109"
clear all
set more off
adopath++ "ado\"
cdate

*************************
******** EXAMPLE ********
*************************
mkpath savings\graphs\example

* Example 1
foreach spec in fe nofe {
	estimates use savings/estimates/example/networth_`spec'_1
	clear
	set obs 18
	gen year = _n+1994
	expand 2, gen(t)
	label var year "Year"

	predict networth
	predict se, stdp

	gen ciup = networth + 1.96*se
	gen cilo = networth - 1.96*se

	twoway  (connected networth year if t==1, color(black)) ///
					(rline ciup cilo year if t==1, lp(-) color(black)) ///
					(connected networth year if t==0, color(gs11)) ///
					(rline ciup cilo year if t==0, lp(-) color(gs11)) ///
					, graphr(color(white)) legend(order(1 "Inheritance in 2000" 3 "Inheritance in 2006")) ///
					ylabel(-1(0.5)2.5) xlabel(1994(3)2012) ytitle("Years of perm. income") ///
					xline(1999.5, lp(-) lc(black)) xline(2005.5, lp(-) lc(gs11)) xtitle("") ysize(3)
	graph export savings\graphs\example\avgs_`spec'_1.pdf, replace
	graph export savings\graphs\example\avgs_`spec'_1.png, replace


	gen difference = 0
	gen ciup_d = .
	gen cilo_d = .
	foreach y of numlist 1995/1998 2000/2012 {
		qui: replace difference = _b[1.t#`y'.year] if year==`y'
		qui: replace ciup_d = _b[1.t#`y'.year] + 1.96*_se[1.t#`y'.year] if year==`y'
		qui: replace cilo_d = _b[1.t#`y'.year] - 1.96*_se[1.t#`y'.year] if year==`y'
	}

	twoway 	(scatter difference year if year<2006 & t, color(black)) ///
					(rcap ciup_d cilo_d year if year<2006 & t, lc(black)) ///
					(scatter difference year if year>=2006 & t, color(gs11)) ///
					(rcap ciup_d cilo_d year if year>=2006 & t, lc(gs11)) ///
					, graphr(color(white)) legend(order(1 "Difference between the series")) ///
					xlabel(1994(3)2012) xline(1999.5, lp(-) lc(black)) ytitle("Years of perm. income") ///
					xline(2005.5, lp(-) lc(gs11)) ///
					ylabel(-1(0.5)1) yline(0, lc(black) lw(medthick)) xtitle("") ysize(3)
	graph export savings\graphs\example\diffs_`spec'_1.pdf, replace 
	graph export savings\graphs\example\diffs_`spec'_1.png, replace 
}

* Example 2
foreach spec in fe nofe {
	estimates use savings/estimates/example/networth_`spec'_2
	clear
	set obs 18
	gen year = _n+1994
	expand 2, gen(t)
	label var year "Year"

	predict networth
	predict se, stdp

	gen ciup = networth + 1.96*se
	gen cilo = networth - 1.96*se

	twoway  (connected networth year if t==1, color(black)) ///
					(rline ciup cilo year if t==1, lp(-) color(black)) ///
					(connected networth year if t==0, color(gs11)) ///
					(rline ciup cilo year if t==0, lp(-) color(gs11)) ///
					, graphr(color(white)) legend(order(1 "Inheritance in 2001" 3 "Inheritance in 2010")) ///
					ylabel(-1(0.5)3) xlabel(1994(3)2012) ytitle("Years of perm. income") ///
					xline(1999.5, lp(-) lc(black)) xline(2009.5, lp(-) lc(gs11)) xtitle("")  ysize(3)
	graph export savings\graphs\example\avgs_`spec'_2.pdf, replace
	graph export savings\graphs\example\avgs_`spec'_2.png, replace

	
	gen difference = 0
	gen ciup_d = .
	gen cilo_d = .
	foreach y of numlist 1995/1998 2000/2012 {
		qui: replace difference = _b[1.t#`y'.year] if year==`y'
		qui: replace ciup_d = _b[1.t#`y'.year] + 1.96*_se[1.t#`y'.year] if year==`y'
		qui: replace cilo_d = _b[1.t#`y'.year] - 1.96*_se[1.t#`y'.year] if year==`y'
	}

	twoway 	(scatter difference year if year<2010 & t, color(black)) ///
					(rcap ciup_d cilo_d year if year<2010 & t, lc(black)) ///
					(scatter difference year if year>=2010 & t, color(gs11)) ///
					(rcap ciup_d cilo_d year if year>=2010 & t, lc(gs11)) ///
					, graphr(color(white)) legend(order(1 "Difference between the series")) ///
					xlabel(1994(3)2012) ytitle("Years of perm. income") ///
					xline(1999.5, lp(-) lc(black)) xline(2009.5, lp(-) lc(gs11)) ///
					ylabel(-1(0.5)1.5) yline(0, lc(black) lw(medthick)) xtitle("") ysize(3)
	graph export savings\graphs\example\diffs_`spec'_2.pdf, replace 
	graph export savings\graphs\example\diffs_`spec'_2.png, replace 
}

* translate into event study
foreach example in 1 2 {
	foreach type in event eventrestr {
		clear
		set obs 15
		gen time = _n + 94
		forval i = 95/109 {
			local real = `i'-100
			label define timef `i' "`real'", add
		}
		label val time timef
		
		estimates use savings/estimates/example/networth_`type'_`example'
		
		gen beta = .
		gen ciup = .
		gen cilo = .
		replace beta = 0 if time ==99
		foreach i of numlist 96/98 100/109 {
				qui: cap replace beta = _b[`i'.timefid] if time == `i'
				qui: cap replace ciup = _b[`i'.timefid] + 1.96*_se[`i'.timefid] if time == `i'
				qui: cap replace cilo = _b[`i'.timefid] - 1.96*_se[`i'.timefid] if time == `i'
		}
			
		twoway (scatter beta time, mc(black)) ///
					 (rcap ciup cilo time, color(black)) ///
					 , graphr(color(white)) xtitle("Years from parental death") ///
					 xlabel(95(1)109, val) ylabel(-0.4(0.2)1.4) ytitle("Years of perm. income") ///
					 xline(100, lp(-) lc(black)) yline(0, lw(medthick) lc(black)) ///
					 legend(off) ysize(3)
		graph export savings\graphs\example/`type'_`example'.pdf, replace
		graph export savings\graphs\example/`type'_`example'.png, replace
	}
}



******************************
******** FN --> EVENT ********
******************************
mkpath savings\graphs\FN_event

foreach var in networth liqassets {
	foreach type in FN {
		foreach control in yearonly_balance yearcohort_balance yearonly_nobalance yearcohort_nobalance {
			clear
			set obs 12
			gen time = _n + 96
			forval i = 97/108 {
				local real = `i'-100
				label define timef `i' "`real'", add
			}
			label val time timef
			
			estimates use savings/estimates/FN/`var'_`type'_`control'
			
			if "`type'" == "FN" {
				local toplot "1.treat#"
			}
			else {
				local toplot ""
			}
			
			gen beta = .
			gen ciup = .
			gen cilo = .
			replace beta = 0 if time ==99
			foreach i of numlist 97/98 100/108 {
					qui: cap replace beta = _b[`toplot'`i'.timefid] if time == `i'
					qui: cap replace ciup = _b[`toplot'`i'.timefid] + 1.96*_se[`toplot'`i'.timefid] if time == `i'
					qui: cap replace cilo = _b[`toplot'`i'.timefid] - 1.96*_se[`toplot'`i'.timefid] if time == `i'
			}
				
			twoway (scatter beta time, mc(black)) ///
						 (rcap ciup cilo time, color(black)) ///
						 , graphr(color(white)) xtitle("Years from parental death") ///
						 xlabel(97(1)108, val) ylabel(-0.2(0.2)1.2) ///
						 xline(100, lp(-) lc(black)) yline(0, lw(medthick) lc(black)) ///
						 legend(off) ysize(3) ytitle("Years of perm. income")
			graph export savings\graphs\FN_event/`var'_`type'_`control'.pdf, replace
			graph export savings\graphs\FN_event/`var'_`type'_`control'.png, replace
		}
	}
}
foreach var in networth liqassets {
	foreach type in event {
		foreach control in yearonly yearcohort {
			clear
			set obs 12
			gen time = _n + 96
			forval i = 97/108 {
				local real = `i'-100
				label define timef `i' "`real'", add
			}
			label val time timef
			
			estimates use savings/estimates/FN/`var'_`type'_`control'
			
			if "`type'" == "FN" {
				local toplot "1.treat#"
			}
			else {
				local toplot ""
			}
			
			gen beta = .
			gen ciup = .
			gen cilo = .
			replace beta = 0 if time ==99
			foreach i of numlist 97/98 100/108 {
					qui: cap replace beta = _b[`toplot'`i'.timefid] if time == `i'
					qui: cap replace ciup = _b[`toplot'`i'.timefid] + 1.96*_se[`toplot'`i'.timefid] if time == `i'
					qui: cap replace cilo = _b[`toplot'`i'.timefid] - 1.96*_se[`toplot'`i'.timefid] if time == `i'
			}
				
			twoway (scatter beta time, mc(black)) ///
						 (rcap ciup cilo time, color(black)) ///
						 , graphr(color(white)) xtitle("Years from parental death") ///
						 xlabel(97(1)108, val) ylabel(-0.2(0.2)1.2) ///
						 xline(100, lp(-) lc(black)) yline(0, lw(medthick) lc(black)) ///
						 legend(off) ysize(3) ytitle("Years of perm. income")
			graph export savings\graphs\FN_event/`var'_`type'_`control'.pdf, replace
			graph export savings\graphs\FN_event/`var'_`type'_`control'.png, replace
		}
	}
}

**********************************
***** REDUCED  FORM - LEVELS *****
**********************************
mkpath savings\graphs\rf_levels
/* Networth */
	** Create dataset **
foreach var in networth liqassets {
	clear all
	*est replay

	set obs 15
	gen time = _n + 94
	forval i = 95/109 {
		local real = `i'-100
		label define timef `i' "`real'", add
	}
	label val time timef
	
	foreach inher in 0 2 {
		
		estimates use savings/estimates/rev1/`var'_`inher'_main
		
		gen beta_`inher' = .
		gen ciup_`inher' = .
		gen cilo_`inher' = .
		foreach i of numlist 95/98 99 100/109 {
			quietly {
				if `i'==99 {
					replace beta_`inher' = 0 if time == `i'
					replace ciup_`inher' = . if time == `i'
					replace cilo_`inher' = . if time == `i'
				}
				else{
					replace beta_`inher' = _b[`i'.timefid] if time == `i'
					replace ciup_`inher' = _b[`i'.timefid] + 1.96*_se[`i'.timefid] if time == `i'
					replace cilo_`inher' = _b[`i'.timefid] - 1.96*_se[`i'.timefid] if time == `i'
				}
			}
		}
	}

	** Plot levels
	foreach inher in 0 2 {
		twoway (scatter beta_`inher' time, mc(black)) ///
					 (rcap ciup_`inher' cilo_`inher' time, color(black)) ///
					 , graphr(color(white)) xtitle("Years from parental death") ///
					 xlabel(95(1)109, val) ylabel(-0.2(0.2)1) ///
					 xline(100, lp(-) lc(black)) yline(0, lw(medthick) lc(black)) ///
					 legend(off) ytitle("Years of perm. income")
		graph export savings\graphs\rf_levels/`var'_`inher'_main.pdf, replace
		graph export savings\graphs\rf_levels/`var'_`inher'_main.png, replace
	}
}	
/* Wealth components, normalized */
	** Create dataset **

	clear all

	set obs 15
	gen time = _n + 94
	forval i = 95/109 {
		local real = `i'-100
		label define timef `i' "`real'", add
	}
	label val time timef
	
	** load and save reg results ** 
	foreach var in liqassets hequity finw liqdebts {
		foreach inher in 0 2 {
		
			estimates use savings/estimates/rev1/`var'_`inher'_main
			
			gen `var'_beta_`inher' = .
			gen `var'_ciup_`inher' = .
			gen `var'_cilo_`inher' = .
			foreach i of numlist 95/98 99 100/109 {
				quietly {
					if `i'==99 {
						replace `var'_beta_`inher' = 0 if time == `i'
						replace `var'_ciup_`inher' = . if time == `i'
						replace `var'_cilo_`inher' = . if time == `i'
					}
					else{
						replace `var'_beta_`inher' = _b[`i'.timefid] if time == `i'
						replace `var'_ciup_`inher' = _b[`i'.timefid] + 1.96*_se[`i'.timefid] if time == `i'
						replace `var'_cilo_`inher' = _b[`i'.timefid] - 1.96*_se[`i'.timefid] if time == `i'
					}
				}
			}
		}
	}
	
	gen base = 0
	** Modify for rarea **
	foreach inher in 0 2 {
		replace liqassets_beta_`inher' = liqassets_beta_`inher' + finw_beta_`inher' + hequity_beta_`inher' - liqdebts_beta_`inher'
		replace finw_beta_`inher' = finw_beta_`inher' + hequity_beta_`inher' - liqdebts_beta_`inher'
		replace hequity_beta_`inher' = hequity_beta_`inher' - liqdebts_beta_`inher'
		replace liqdebts_beta_`inher' = - liqdebts_beta_`inher'
	}

	** Plot ratios **
	
	foreach inher in 0 2 {
		twoway (rarea liqassets_beta_`inher' base time, color(gs6)) ///
					 (rarea finw_beta_`inher' base time, color(gs10)) ///
					 (rarea hequity_beta_`inher' base time, color(gs14)) ///
					 (rarea liqdebts_beta_`inher' base time, color(gs2)) ///
					 if time>98, graphr(color(white))  xtitle("Years from parental death") ///
					 legend(region(color(white)) order( 1 "Liquid assets" ///
					 2 "Financial assets" 3 "Housing equity" ///
					 4 "Non-collateralized debts (reversed)")) ///
					xlabel(99(1)109, val) ylabel(0(0.2)1) yscale(range(-0.1 1)) ytitle("Years of perm. income")
		graph export savings\graphs\rf_levels\wealthcomp_`inher'_fe_main.pdf, replace
		graph export savings\graphs\rf_levels\wealthcomp_`inher'_fe_main.png, replace width(1500)
	}

*
*************************************
************ APPENDIX ***************
*************************************

*** - BY INHERITANCE SIZE ***

clear all

set obs 15
gen time = _n + 94
forval i = 95/109 {
	local real = `i'-100
	label define timef `i' "`real'", add
}
label val time timef

foreach inher in 1 2 3 {
	
	estimates use savings/estimates/rev1/networth_2_inhgr`inher'
	
	gen beta_`inher' = .
	gen ciup_`inher' = .
	gen cilo_`inher' = .
	foreach i of numlist 95/98 99 100/109 {
		quietly {
			if `i'==99 {
				replace beta_`inher' = 0 if time == `i'
				replace ciup_`inher' = . if time == `i'
				replace cilo_`inher' = . if time == `i'
			}
			else{
				replace beta_`inher' = _b[`i'.timefid] if time == `i'
				replace ciup_`inher' = _b[`i'.timefid] + 1.96*_se[`i'.timefid] if time == `i'
				replace cilo_`inher' = _b[`i'.timefid] - 1.96*_se[`i'.timefid] if time == `i'
			}
		}
	}
	
	gen beta_`inher'_norm = beta_`inher'/_b[101.timefid]
	gen ciup_`inher'_norm = ciup_`inher'/_b[101.timefid]
	gen cilo_`inher'_norm = cilo_`inher'/_b[101.timefid]
	
}
*
cap drop time2 time3
gen time2 = time - 0.15
gen time3 = time + 0.15
twoway 	(connected beta_2 time, color(gs8) lp(.) ms(T)) ///
		(rcap ciup_2 cilo_2 time, color(gs8) lp(.)) ///
		(connected beta_1 time2, color(gs12) lp(-) ms(O)) ///
		(rcap ciup_1 cilo_1 time2, color(gs12) lp(-)) ///
		(connected beta_3 time3, color(black) ms(Oh)) ///
		(rcap ciup_3 cilo_3 time3, color(black) ) ///
	, graphr(color(white)) ///
	xlabel(95(1)109, val)  ytitle("Years of perm. income") ///
	legend(order(- "Imp. inheritance:" 3 "1 to 2.5" 1 "2.5 to 5" 5 "> 5") ///
		region(color(white)) r(1)) xtitle("Years from parental death") ///
	xline(100, lp(-) lc(black)) yline(0, lw(medthick) lc(black))
graph export savings\graphs\rf_levels\robustness_byinher.pdf, replace
	

twoway 	(connected beta_2_norm time, color(gs8) lp(.) ms(T)) ///
		(connected beta_1_norm time, color(gs12) lp(-) ms(O)) ///
		(connected beta_3_norm time, color(black) ms(Oh)) ///
	, graphr(color(white)) ///
	xlabel(95(1)109, val)  ytitle("Wealth change ratio") ///
	legend(order(- "Imp. inheritance:" 2 "1 to 2.5" 1 "2.5 to 5" 3 "> 5") ///
		region(color(white)) r(1)) xtitle("Years from parental death") ///
	xline(100, lp(-) lc(black)) yline(0, lw(medthick) lc(black)) 
graph export savings\graphs\rf_levels\robustness_byinher_norm.pdf, replace


*** COMPARE ABSOLUTE AND NORMALIZED EVOLUTION ***
clear all

set obs 15
gen time = _n + 94
forval i = 95/109 {
	local real = `i'-100
	label define timef `i' "`real'", add
}
label val time timef

foreach var in networth networth_abs {
	
	estimates use savings/estimates/rev1/`var'_2_main
	
	gen beta_`var' = .
	gen ciup_`var' = .
	gen cilo_`var' = .
	
	foreach i of numlist 95/98 99 100/109 {
		quietly {
			if `i'==99 {
				replace beta_`var' = 0 if time == `i'
				replace ciup_`var' = . if time == `i'
				replace cilo_`var' = . if time == `i'
			}
			else{
				replace beta_`var' = _b[`i'.timefid] if time == `i'
				replace ciup_`var' = _b[`i'.timefid] + 1.96*_se[`i'.timefid] if time == `i'
				replace cilo_`var' = _b[`i'.timefid] - 1.96*_se[`i'.timefid] if time == `i'
			}
		}
	}
	
	gen beta_`var'_norm = beta_`var'/_b[101.timefid]
	gen ciup_`var'_norm = ciup_`var'/_b[101.timefid]
	gen cilo_`var'_norm = cilo_`var'/_b[101.timefid]
	
}

cap drop time2 time3
gen time2 = time + 0.15
gen time3 = time - 0.15

twoway	(scatter beta_networth_abs_norm time, m(i)) ///
		(scatter beta_networth_abs_norm time2, color(gs8) ms(Oh)) ///
		(rcap ciup_networth_abs_norm cilo_networth_abs_norm time2, color(gs8) lp(-)) ///
		(scatter beta_networth_norm time3, color(black) ms(O)) ///
		(rcap ciup_networth_norm cilo_networth_norm time3, color(black) ) ///
	, graphr(color(white)) ysize(5) xsize(8) ///
	xlabel(95(1)109, val)  ytitle("") ///
	legend(order(4 "Normalized by permanent income" 2 "Absolute values" ) ///
		region(color(white)) r(1)) xtitle("Years from parental death") ///
	xline(100, lp(-) lc(black)) yline(0, lw(medthick) lc(black))
graph export savings\graphs\rf_levels\robustness_norm_VS_abs.pdf, replace


*** COMPARE ACTUAL MEASURE OF WEALTH WITH ONE WHERE HOUSE VALUES ARE INFLATED 20% ***
clear all

set obs 15
gen time = _n + 94
forval i = 95/109 {
	local real = `i'-100
	label define timef `i' "`real'", add
}
label val time timef

foreach var in networth networth_infl {
	
	estimates use savings/estimates/rev1/`var'_2_main
	
	gen beta_`var' = .
	gen ciup_`var' = .
	gen cilo_`var' = .
	
	foreach i of numlist 95/98 99 100/109 {
		quietly {
			if `i'==99 {
				replace beta_`var' = 0 if time == `i'
				replace ciup_`var' = . if time == `i'
				replace cilo_`var' = . if time == `i'
			}
			else{
				replace beta_`var' = _b[`i'.timefid] if time == `i'
				replace ciup_`var' = _b[`i'.timefid] + 1.96*_se[`i'.timefid] if time == `i'
				replace cilo_`var' = _b[`i'.timefid] - 1.96*_se[`i'.timefid] if time == `i'
			}
		}
	}
	
	gen beta_`var'_norm = beta_`var'/_b[101.timefid]
	gen ciup_`var'_norm = ciup_`var'/_b[101.timefid]
	gen cilo_`var'_norm = cilo_`var'/_b[101.timefid]
	
}


cap drop time2 
cap drop time3
gen time2 = time + 0.15
gen time3 = time - 0.15

twoway	(scatter beta_networth_infl_norm time, m(i)) ///
		(scatter beta_networth_norm time3, color(black) ms(O)) ///
		(rcap ciup_networth_norm cilo_networth_norm time3, color(black) ) ///
		(scatter beta_networth_infl_norm time2, color(gs8) ms(Oh)) ///
		(rcap ciup_networth_infl_norm cilo_networth_infl_norm time2, color(gs8) lp(-)) ///
	, graphr(color(white)) ysize(5) xsize(8) ///
	xlabel(95(1)109, val)  ytitle("") ///
	legend(order(2 "Baseline" 4 "Housing values inflated by 20%") ///
		region(color(white)) r(1)) xtitle("Years from parental death") ///
	xline(100, lp(-) lc(black)) yline(0, lw(medthick) lc(black))
graph export savings\graphs\rf_levels\robustness_base_VS_inflated_housing.pdf, replace


*** COMPARE MAIN RESULTS AT THE INDIVIDUAL AND HOUSEHOLD LEVEL ***
clear all

set obs 15
gen time = _n + 94
forval i = 95/109 {
	local real = `i'-100
	label define timef `i' "`real'", add
}
label val time timef

foreach var in networth_2_hh_sample hh_networth_2_main liqworth_2_hh_sample hh_liqworth_2_main {
	
	estimates use savings/estimates/rev1/`var'
	
	gen beta_`var' = .
	gen ciup_`var' = .
	gen cilo_`var' = .
	
	foreach i of numlist 95/98 99 100/109 {
		quietly {
			if `i'==99 {
				replace beta_`var' = 0 if time == `i'
				replace ciup_`var' = . if time == `i'
				replace cilo_`var' = . if time == `i'
			}
			else{
				replace beta_`var' = _b[`i'.timefid] if time == `i'
				replace ciup_`var' = _b[`i'.timefid] + 1.96*_se[`i'.timefid] if time == `i'
				replace cilo_`var' = _b[`i'.timefid] - 1.96*_se[`i'.timefid] if time == `i'
			}
		}
	}
	
	gen beta_`var'_norm = beta_`var'/_b[101.timefid]
	gen ciup_`var'_norm = ciup_`var'/_b[101.timefid]
	gen cilo_`var'_norm = cilo_`var'/_b[101.timefid]
	
}


cap drop time2 
cap drop time3
gen time2 = time + 0.15
gen time3 = time - 0.15

twoway	(scatter beta_hh_networth_2_main_norm time, m(i)) ///
		(scatter beta_networth_2_hh_sample_norm time3, color(black) ms(O)) ///
		(rcap ciup_networth_2_hh_sample_norm cilo_networth_2_hh_sample_norm time3, color(black) ) ///
		(scatter beta_hh_networth_2_main_norm time2, color(gs8) ms(Oh)) ///
		(rcap ciup_hh_networth_2_main_norm cilo_hh_networth_2_main_norm time2, color(gs8) lp(-)) ///
	, graphr(color(white)) ysize(5) xsize(8) ///
	xlabel(95(1)109, val)  ytitle("") ///
	legend(order(2 "Individual level" 4 "Household level") ///
		region(color(white)) r(1)) xtitle("Years from parental death") ///
	xline(100, lp(-) lc(black)) yline(0, lw(medthick) lc(black))
	graph export savings\graphs\rf_levels\robustness_networth_individual_VS_hhd.pdf, replace
	
twoway	(scatter beta_hh_liqworth_2_main_norm time, m(i)) ///
		(scatter beta_liqworth_2_hh_sample_norm time3, color(black) ms(O)) ///
		(rcap ciup_liqworth_2_hh_sample_norm cilo_liqworth_2_hh_sample_norm time3, color(black) ) ///
		(scatter beta_hh_liqworth_2_main_norm time2, color(gs8) ms(Oh)) ///
		(rcap ciup_hh_liqworth_2_main_norm cilo_hh_liqworth_2_main_norm time2, color(gs8) lp(-)) ///
	, graphr(color(white)) ysize(5) xsize(8) ///
	xlabel(95(1)109, val)  ytitle("") ///
	legend(order(2 "Individual level" 4 "Household level") ///
		region(color(white)) r(1)) xtitle("Years from parental death") ///
	xline(100, lp(-) lc(black)) yline(0, lw(medthick) lc(black))
	graph export savings\graphs\rf_levels\robustness_liqworth_individual_VS_hhd.pdf, replace


*** COMPARE MAIN RESULTS FOR LIQUID AND ILLIQUID INHERITANCES ***
clear all

set obs 15
gen time = _n + 94
forval i = 95/109 {
	local real = `i'-100
	label define timef `i' "`real'", add
}
label val time timef

foreach var in networth_2_sire_liquid networth_2_sire_illiquid {
	
	estimates use savings/estimates/rev1/`var'
	
	gen beta_`var' = .
	gen ciup_`var' = .
	gen cilo_`var' = .
	
	foreach i of numlist 95/98 99 100/109 {
		quietly {
			if `i'==99 {
				replace beta_`var' = 0 if time == `i'
				replace ciup_`var' = . if time == `i'
				replace cilo_`var' = . if time == `i'
			}
			else{
				replace beta_`var' = _b[`i'.timefid] if time == `i'
				replace ciup_`var' = _b[`i'.timefid] + 1.96*_se[`i'.timefid] if time == `i'
				replace cilo_`var' = _b[`i'.timefid] - 1.96*_se[`i'.timefid] if time == `i'
			}
		}
	}
	
	gen beta_`var'_n = beta_`var'/_b[101.timefid]
	gen ciup_`var'_n = ciup_`var'/_b[101.timefid]
	gen cilo_`var'_n = cilo_`var'/_b[101.timefid]
	
}


cap drop time2 
cap drop time3
gen time2 = time + 0.15
gen time3 = time - 0.15

twoway	///
		(scatter beta_networth_2_sire_liquid_n time3, color(black) ms(O)) ///
		(rcap ciup_networth_2_sire_liquid_n cilo_networth_2_sire_liquid_n time3, color(black) ) ///
		(scatter beta_networth_2_sire_illiquid_n time2, color(gs8) ms(Oh)) ///
		(rcap ciup_networth_2_sire_illiquid_n cilo_networth_2_sire_illiquid_n time2, color(gs8) lp(-)) ///
	, graphr(color(white)) ysize(5) xsize(8) ///
	xlabel(95(1)109, val)  ytitle("") ///
	legend(order(2 "Liquid inheritance" 4 "Illiquid inheritance") ///
		region(color(white)) r(1)) xtitle("Years from parental death") ///
	xline(100, lp(-) lc(black)) yline(0, lw(medthick) lc(black))
	graph export savings\graphs\rf_levels\robustness_networth_liquid_vs_illiquid.pdf, replace
	
	
	
**********************************
*****   FOR  PRESENTATIONS   *****
***** REDUCED  FORM - LEVELS *****
**********************************
mkpath savings\graphs\presentations
/* Networth */
	** Create dataset **
foreach var in networth liqassets hequity finw n_owned {
	foreach sample in main young exp_young old exp_old costr exp_costr uncostr exp_uncostr {
		clear all
		*est replay

		set obs 16
		gen time = _n + 93
		label define timef 94 "<-5"
		forval i = 95/109 {
			local real = `i'-100
			label define timef `i' "`real'", add
		}
		label val time timef
		
		foreach inher in 2 {
			
			estimates use savings/estimates/levels_reducedform/`var'_`inher'_fe_`sample'
			
			gen beta_`inher' = _b[0.inhgroup_all]
			gen ciup_`inher' = _b[0.inhgroup_all] + 1.96*_se[0.inhgroup_all]
			gen cilo_`inher' = _b[0.inhgroup_all] - 1.96*_se[0.inhgroup_all]
			foreach i of numlist 95/98 99 100/109 {
				quietly {
					if `i'==99 {
						replace beta_`inher' = 0 if time == `i'
						replace ciup_`inher' = . if time == `i'
						replace cilo_`inher' = . if time == `i'
					}
					else{
						replace beta_`inher' = _b[`i'.inhgroup_all] if time == `i'
						replace ciup_`inher' = _b[`i'.inhgroup_all] + 1.96*_se[`i'.inhgroup_all] if time == `i'
						replace cilo_`inher' = _b[`i'.inhgroup_all] - 1.96*_se[`i'.inhgroup_all] if time == `i'
					}
				}
			}
		}

		** Plot levels
		foreach inher in 2 {
			twoway (scatter beta_`inher' time, mc(black)) ///
						 (rcap ciup_`inher' cilo_`inher' time, color(black)) ///
						 , graphr(color(white)) xtitle("Years from parental death") ///
						 xlabel(94(1)109, val) ylabel(-0.4(0.2)1.2) ///
						 xline(100, lp(-) lc(black)) yline(0, lw(medthick) lc(black)) ///
						 legend(off)
			graph export savings\graphs\presentations/`var'_`inher'_fe_`sample'.pdf, replace
			graph export savings\graphs\presentations/`var'_`inher'_fe_`sample'.png, replace
		}
	}
}	
*
	** Create dataset **
foreach var in networth liqassets hequity finw n_owned {
	foreach sample in main {
		clear all
		*est replay

		set obs 16
		gen time = _n + 93
		label define timef 94 "<-5"
		forval i = 95/109 {
			local real = `i'-100
			label define timef `i' "`real'", add
		}
		label val time timef
		
		foreach inher in 0 {
			
			estimates use savings/estimates/levels_reducedform/`var'_`inher'_fe_`sample'
			
			gen beta_`inher' = _b[0.inhgroup_all]
			gen ciup_`inher' = _b[0.inhgroup_all] + 1.96*_se[0.inhgroup_all]
			gen cilo_`inher' = _b[0.inhgroup_all] - 1.96*_se[0.inhgroup_all]
			foreach i of numlist 95/98 99 100/109 {
				quietly {
					if `i'==99 {
						replace beta_`inher' = 0 if time == `i'
						replace ciup_`inher' = . if time == `i'
						replace cilo_`inher' = . if time == `i'
					}
					else{
						replace beta_`inher' = _b[`i'.inhgroup_all] if time == `i'
						replace ciup_`inher' = _b[`i'.inhgroup_all] + 1.96*_se[`i'.inhgroup_all] if time == `i'
						replace cilo_`inher' = _b[`i'.inhgroup_all] - 1.96*_se[`i'.inhgroup_all] if time == `i'
					}
				}
			}
		}

		** Plot levels
		foreach inher in 0 {
			twoway (scatter beta_`inher' time, mc(black)) ///
						 (rcap ciup_`inher' cilo_`inher' time, color(black)) ///
						 , graphr(color(white)) xtitle("Years from parental death") ///
						 xlabel(94(1)109, val) ylabel(-0.4(0.2)1.2) ///
						 xline(100, lp(-) lc(black)) yline(0, lw(medthick) lc(black)) ///
						 legend(off)
			graph export savings\graphs\presentations/`var'_`inher'_fe_`sample'.pdf, replace
			graph export savings\graphs\presentations/`var'_`inher'_fe_`sample'.png, replace
		}
	}
}	
*
