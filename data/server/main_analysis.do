cd "D:\Data\workdata\705109"
clear all
set more off
adopath++ "ado\"
cdate

timer clear 1
timer on 1

*****************************
********* LOAD DATA *********
*****************************
use data/workdata, clear
drop if inheritance==.
*keep if minobs<=-2 & maxobs>=2

gen ln_dispinc = log(max(1,dispinc))

/* Measures at shock */
bys id (year): egen alder_atshock = total(alder*(timefromshock==0))
bys id (year): egen sire_alder_atshock = total(sire_alder*(timefromshock==-1))
	replace sire_alder_atshock = sire_alder_atshock +1
	replace sire_alder_atshock =. if sire_alder_atshock==1
drop if alder_atshock==0
keep if alder_atshock<50

/* is there in timefromschok-1? */
by id: egen x = max(timefromshock==-1)
drop if x==0

**********************************************************
**** Compute variance of shocks for model calibration ****
**********************************************************
	
cap drop ln_dispinc
cap drop incvar 
gen incvar = dispinc
gen ln_dispinc = log(max(1,incvar)) if dispinc!=.
xtset id year

cap drop lninc d_lninc*
cap drop sum_d_*
reg ln_dispinc b40.alder if incvar>10000
predict lninc, r, if incvar>10000

gen d_lninc = d.lninc
gen d_lninc_m1 = l.d_lninc
gen d_lninc_p1 = f.d_lninc

gen sum_d_ln = d_lninc + d_lninc_m1 + d_lninc_p1

corr d_lninc d_lninc_p1 , cov
di "Sigma_epsilon: " `=sqrt(-r(cov_12))'

corr d_lninc sum_d_ln , cov
di "Sigma_phi: " `=sqrt(r(cov_12))'

** Write moments to file **
file open toexport using toexport/income_stds.csv, write
file write toexport "sigma_epsilon, " (`=sqrt(-r(cov_12))') _n "sigma_phi, " (`=sqrt(r(cov_12))')
file close toexport

*****************************
******** SET GLOBALS ********
*****************************
global outcomes networth liqworth illiqworth liqassets liqdebts finw hequity hvalue mortgage hequity_infl networth_infl
global outcomes_abs	networth_abs liqworth_abs illiqworth_abs liqassets_abs liqdebts_abs finw_abs hequity_abs hvalue_abs mortgage_abs
global other_outcomes n_owned home_own multi_own married antboernh dispinc labor_income salary has_labor has_qpripen qpripen_abs qpripen_norm has_qarbpen qarbpen_abs qarbpen_norm
global spo_outcomes spo_dispinc spo_has_labor spo_salary spo_labor_income spo_networth spo_n_owned 
global household_outcomes hh_dispinc hh_labor_income hh_liqworth hh_networth
global mainfes i.year#i.cohort


******************************
*** CREATE ADDITIONAL VARS ***
******************************

/* Asset grouping */ 
gen liqworth = liqassets - liqdebts
gen illiqworth = hequity + finw 


/* Inflated net worth variable - robustness check */
gen hequity_infl = hequity + hvalue*0.2 
gen networth_infl = networth + hvalue*0.2


/* Get back to base values*/
foreach var in $outcomes {
	gen `var'_abs = `var'*perminc/1000
}
replace dispinc = dispinc/1000
replace labor_income = labor_income/1000
replace salary = salary/1000

gen inhgroup = inheritance > 1/12
replace inhgroup = 2 if inheritance>1 			// sample of interest (==2)
gen inheritance_cens = min(inheritance, 15)		// for robustness check

replace timefromshock = timefromshock + 100
cap drop identifiable
gen identifiable = timefromshock>=95			// Double normalization necessary for identification

gen timefid = timefromshock*identifiable

cap drop inhgroup_all
gen inhgroup_all = timefromshock*identifiable
forval i = 0/2 {
	cap drop inhgroup_`i'
	gen inhgroup_`i' = timefromshock*identifiable*(inhgroup==`i')
}
*
bys id: egen lq_constr = max((liqassets<1/12)*(timefromshock==99))

cap drop home_own
gen home_own = hvalue>0
cap drop owner_before
bys id: egen owner_before = total(inrange(timefromshock,97,99)*home_own)
gen multi_own = n_owned>1 & n_owned!=.

gen _x = spo_pnr if timefid==99
gen spo_pnr_atshock = ""
bys id (_x): replace spo_pnr_atshock = _x[_N]
gen consistent_hh = spo_pnr_atshock == spo_pnr
drop _x


replace spo_perminc = 0 if spo_pnr ==""
replace spo_networth = spo_networth/spo_perminc
replace spo_liqdebts = spo_liqdebts/spo_perminc
replace spo_liqassets = spo_liqassets/spo_perminc
replace spo_hequity = spo_hequity/spo_perminc

replace spo_dispinc = max(0, spo_dispinc)/1000
replace spo_labor_income = max(0, spo_labor_income)/1000
replace spo_salary = max(0, spo_salary)/1000
gen spo_has_labor = spo_labor_income>0

foreach var in $spo_outcomes {
	replace `var' = . if married !=1 | consistent_hh==0
}


cap drop x
gen x = 100*sire_hvalue/(sire_totassets)
bys id: egen sire_illiquid = max(x*(timefid==99))

foreach var in networth liqdebts liqassets {
	replace spo_`var' = 0 if spo_pnr ==""
	cap drop hh_`var'
	gen hh_`var' = (`var'*perminc + spo_`var'*spo_perminc)/(perminc + spo_perminc)
}
cap drop hh_liqworth
gen hh_liqworth = hh_liqassets - hh_liqdebts
foreach var in dispinc labor_income {
	replace spo_`var' = 0 if spo_pnr ==""
	cap drop hh_`var'
	gen hh_`var' = (`var' + spo_`var')
}
gen cohort = year(foed_dag)

cap drop married
gen married = spo_pnr!=""

replace labor_income = max(0, labor_income)
replace salary = max(0, salary)
gen has_labor = labor_income>0

foreach pens in qpripen qarbpen {
	replace `pens' = max(0, `pens')
	gen has_`pens' = `pens'>0
	gen `pens'_norm = `pens'/perminc
	rename `pens' `pens'_abs
}

fvset base 99 timefromshock inhgroup_all timefid
fvset base 0 inhgroup


**************************************************
******** REDUCED FORM - LEVELS/CATEGORIES ********
**************************************************
// We do not estimate times after 9yrs, too imprecise
keep if inrange(timefromshock,0,+9+100)

// odd quirk of stata 14, recoded year as float and would not accept i. shortcut in $mainfes
// This snippet ensures year is recasted as integer
cap drop x
gen x = year
drop year
rename x year
	
*************************************
/********* Estimation step *********/
*************************************
	/* FULL SAMPLE - main*/
	mkpath savings/estimates/rev1
		foreach var in $spo_outcomes hequity_infl networth_infl $outcomes $outcomes_abs $other_outcomes {
			foreach group in 0 2 {
			reghdfe `var' b99.timefid if inhgroup==`group' & sudden==1, ///
				abs(i.id $mainfes) vce(cluster id)
			estimates save savings/estimates/rev1/`var'_`group'_main, replace
			}
		}
	
	mkpath savings/estimates/rev1
		foreach var in $outcomes $outcomes_abs {
			foreach group in 0 2 {
			reghdfe `var' b99.timefid if inhgroup==`group' & sudden==0, ///
				abs(i.id $mainfes) vce(cluster id)
			estimates save savings/estimates/rev1/`var'_`group'_exp, replace
			}
		}
	
	/* FULL SAMPLE - household outcomes*/
	mkpath savings/estimates/rev1
		foreach var in $household_outcomes {
			foreach group in 2 {
			reghdfe `var' b99.timefid if inhgroup==`group' & sudden==1 & consistent_hh==1, ///
				abs(i.id $mainfes) vce(cluster id)
			estimates save savings/estimates/rev1/`var'_`group'_main, replace
			}
		}		
		cap drop hh_sample
		gen hh_sample = e(sample)
		reghdfe networth b99.timefid if hh_sample, ///
			abs(i.id $mainfes) vce(cluster id)
		estimates save savings/estimates/rev1/networth_2_hh_sample, replace
		reghdfe liqworth b99.timefid if hh_sample, ///
			abs(i.id $mainfes) vce(cluster id)
		estimates save savings/estimates/rev1/liqworth_2_hh_sample, replace	
		
		/* BY LQ COSTRAINTS */
		foreach var in $outcomes {
			foreach group in 2 {
				reghdfe `var' b99.timefid ///
					if inhgroup==`group' & sudden==1 & lq_constr==0, ///
					abs(i.id $mainfes) vce(cluster id)
				estimates save savings/estimates/rev1/`var'_`group'_uncostr, replace
				
				reghdfe `var' b99.timefid ///
					if inhgroup==`group' & sudden==1 & lq_constr==1, ///
					abs(i.id $mainfes) vce(cluster id)
				estimates save savings/estimates/rev1/`var'_`group'_costr, replace
			}
		}
		/* BY INHERITANCE SIZE */
	mkpath savings/estimates/rev1
		foreach var in networth liqworth {
			foreach group in 2 {
			reghdfe `var' b99.timefid if inrange(inheritance_cens,1,2.5) & sudden==1, ///
				abs(i.id $mainfes) vce(cluster id)
			estimates save savings/estimates/rev1/`var'_`group'_inhgr1, replace
			
			reghdfe `var' b99.timefid if inrange(inheritance_cens,2.5,5) & sudden==1, ///
				abs(i.id $mainfes) vce(cluster id)
			estimates save savings/estimates/rev1/`var'_`group'_inhgr2, replace
			
			reghdfe `var' b99.timefid if inrange(inheritance_cens,5,15) & sudden==1, ///
				abs(i.id $mainfes) vce(cluster id)
			estimates save savings/estimates/rev1/`var'_`group'_inhgr3, replace
			}
		}
		/* BY PARENT'S LIQUIDITY */
	mkpath savings/estimates/rev1
		foreach var in networth liqworth {
			foreach group in 2 {
			reghdfe `var' b99.timefid if inhgroup==2 & inrange(sire_illiquid,0,50) & sudden==1, ///
				abs(i.id $mainfes) vce(cluster id)
			estimates save savings/estimates/rev1/`var'_`group'_sire_liquid, replace
	
			reghdfe `var' b99.timefid if inhgroup==2 & inrange(sire_illiquid,50,100) & sudden==1, ///
				abs(i.id $mainfes) vce(cluster id)
			estimates save savings/estimates/rev1/`var'_`group'_sire_illiquid, replace
			}
		}


** TIME RECORDING **
timer off 1
timer list 1
timerec, save(savings/timerec/main_formodels.txt) timer(1)
