cd "D:\Data\workdata\705109"
clear all
set more off
adopath++ "ado\"
cdate

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
*keep if alder_atshock<50

/* is there in timefromschok-1? */
by id: egen x = max(timefromshock==-1)
drop if x==0

** ONLY LARGE INHERITANCES
keep if inheritance>1
keep if sudden ==1

******************************
*** CREATE ADDITIONAL VARS ***
******************************

gen liqworth = liqassets - liqdebts
gen illiqworth = hequity + finw 
gen networth_abs = networth*perminc

replace timefromshock = timefromshock + 100

cap drop identifiable
gen identifiable = timefromshock>=95

gen timefid = timefromshock*identifiable


gen home_own = n_owned>0 if n_owned!=.
gen cohort = year(foed_dag)
gen home_owner = hvalue>0
gen pos_equity = hequity >0
gen pos_liqworth = liqworth>0
gen pos_illiqworth = illiqworth>0

cap drop not_constr
gen not_constr = networth>-0.25


cap drop _sin_1reg
***********************************
*** Compute moments of interest ***
***********************************

cap drop _* 
local varstoexport networth liqworth illiqworth inheritance perminc

foreach var in `varstoexport' {
	cap drop _`var'_mean
	gen _`var'_mean = 0
	
	cap drop x
	cap drop r
	gen x = ln(`var')
	qui: _regress x i.alder b2005.year if `var'>0
	predict r if e(sample), r
	forval a = 25/59 {
		qui: replace _`var'_mean = exp(r)*exp(_b[_cons])*exp(_b[`a'.alder]) if alder==`a' & `var'>0
	}
	cap drop x
	cap drop r
	count if `var'<0
	if r(N)!=0{
		gen x = ln(-`var')
		qui: _regress x i.alder b2005.year if `var'<0
		predict r if e(sample), r
		forval a = 25/59 {
			qui: replace _`var'_mean = -exp(r)*exp(_b[_cons])*exp(_b[`a'.alder]) if alder==`a' & `var'<0
		}
	cap drop x
	cap drop r
	}
}
scalar avg_empi_inh = .8787048
sum _inheritance_mean
scalar avg_impu_inh = r(mean)

replace _inheritance_mean = (_inheritance_mean*avg_empi_inh)/avg_impu_inh
foreach var in `varstoexport' {
	foreach p in 25 50 75 {
		bys alder: egen _`var'_p`p' = DSTpctile(_`var'_mean), p(`p') 
	}
}
foreach p in 25 50 75 {
		bys alder: egen _illiqworth_cond_p`p' = DSTpctile(_illiqworth_mean), p(`p'), if _illiqworth_mean>0
}
gen _illiqworth_cond_mean = _illiqworth_mean if _illiqworth_mean>0
* define scaled vars
gen _pos_liqworth = _liqworth_mean>0
gen _pos_illiqworth = _illiqworth_mean>0

cap drop _not_constr
gen _not_constr = _networth_mean>-0.25




preserve 
	cap drop _constr
	gen _constr = _networth_mean < -0.25
	cap drop _borrow
	gen _borrow = _liqworth_mean<0
	cap drop _neg_netw 
	gen _neg_netw = _networth_mean <0


	gen _corr_norm = .
	gen _corr_abs = .
	levelsof alder, local(alderlevs)
	foreach a of local alderlevs {
		quietly {
			corr perminc networth if alder==`a'
			replace _corr_norm = r(rho) if alder==`a'
			corr perminc networth_abs if alder==`a'
			replace _corr_abs = r(rho) if alder==`a'
		}
	}
	*
	gen count = 1
	collapse (mean) _* (count) count, by(alder)

	order *, alpha
	order alder count

	export delimited using savingstoexport/all.csv, replace
restore

####################
# MOMENTS AT DEATH #
####################

cd "D:\Data\workdata\705109"
clear all
set more off
adopath++ "ado\"
cdate

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
*keep if alder_atshock<50

/* is there in timefromschok-1? */
by id: egen x = max(timefromshock==-1)
drop if x==0

** ONLY LARGE INHERITANCES
keep if inheritance>1

******************************
*** CREATE ADDITIONAL VARS ***
******************************

gen liqworth = liqassets - liqdebts
gen illiqworth = hequity + finw 

/* Temporary - get back to base values*/
foreach var in inheritance networth hequity finw liqassets liqdebts liqworth illiqworth {
	gen `var'_abs = `var'*perminc/1000
}
replace dispinc = dispinc/1000

/* Spouse variables */
egen spouse_id = group(spo_pnr)
bys id: egen spouse_id_base = max(spouse_id*(timefromshock==-1))
foreach var in networth {
	replace spo_`var'=. if spouse_id!=spouse_id_base
	replace spo_`var'=spo_`var'/spo_perminc
	gen spo_`var'_abs = spo_`var'*spo_perminc/1000
}
	replace spo_n_owned=. if spouse_id!=spouse_id_base
*

gen inheritance_cens = min(inheritance,15)
/* Inheritance grouping */
gen 		inhgroup = inheritance>1/12
replace inhgroup = 2 if inheritance>1
*replace inhgroup = 3 if inheritance>2


replace timefromshock = timefromshock + 100

cap drop identifiable
gen identifiable = timefromshock>=95

gen timefid = timefromshock*identifiable

cap drop inhgroup_all
gen inhgroup_all = timefromshock*identifiable
forval i = 0/2 {
	cap drop inhgroup_`i'
	gen inhgroup_`i' = timefromshock*identifiable*(inhgroup==`i')
}
*
bys id: egen lq_constr = max((liqassets<1/12)*(timefromshock==99))
gen old = alder_atshock>40 if alder_atshock!=.
gen married = spouse_id!=0
gen home_own = n_owned>0 if n_owned!=.

gen cohort = year(foed_dag)
gen home_owner = n_own>0 



fvset base 99 timefromshock inhgroup_all timefid
fvset base 0 inhgroup



************************
*** Plot measures at death
***   1) Distributions of age of receipt
***   2) Distribution of age differences
***   3) AVG inheritance over age of receipt
***   4) Distribution inheritance over year of receipt
************************
keep if timefrom==99
replace sire_alder_atshock = sire_alder_atshock+1
gen agediff = sire_alder_atshock - alder_atshock

gen count = 1

hist alder_atshock, disc color(gs6) fintensity(40) ///
	graphr(color(white)) xtitle("Age at parental death")
	graph export savings/moments_jeppe/atdeath/graphs/distributions/ageatdeath.pdf, replace
	preserve 
		collapse (count) count, by(alder_atshock)
		egen total = total(count)
		replace count = count/total
		drop total
		export delimited using savings/toexport/ageatdeath.csv, replace
	restore
	
hist sire_alder_atshock if inrange(sire_alder_atshock,50,95), freq disc color(gs6) fintensity(40) ///
	graphr(color(white)) xtitle("Parent's age at death")
	graph export savings/moments_jeppe/atdeath/graphs/distributions/sire_ageatdeath.pdf, replace
	preserve 
		collapse (count) count if inrange(sire_alder_atshock,50,95), by(sire_alder_atshock)
		egen total = total(count)
		replace count = count/total
		drop total
		export delimited using savings/toexport/sire_ageatdeath.csv, replace
	restore
	
hist agediff if inrange(agediff,18,54), freq disc color(gs6) fintensity(40) ///
	graphr(color(white)) xtitle("Age difference")
	graph export savings/moments_jeppe/atdeath/graphs/distributions/agediff.pdf, replace
	preserve 
		collapse (count) count if inrange(agediff,18,54), by(agediff)
		egen total = total(count)
		replace count = count/total
		drop total
		export delimited using savings/toexport/agediff.csv, replace
	restore
	
	
	foreach p in 10 25 50 75 90 {
		cap drop _pinh_`p'
		bys death_year: egen _pinh_`p' = pctile(inheritance), p(`p'), if inheritance <6
	}
	cap drop pickone
	egen pickone = tag(death_year) if inheritance <6
	
	twoway  (rbar _pinh_50 _pinh_75 death_year ///
					, barw(0.4) bc(none) blc(black) blw(medium) ) ///
					(rbar _pinh_25 _pinh_50 death_year ///
					, barw(0.4) bc(none) blc(black) blw(medium) ) ///
					(rspike _pinh_10 _pinh_25 death_year ///
					, blc(black) blw(medium) ) ///
					(rspike _pinh_90 _pinh_75 death_year ///
					, blc(black) blw(medium) ) ///
					if pickone==1 & death_year<=2010, graphr(color(white)) legend(off)
					
	


/* Correlations */
	// Age difference and age at shock
	binscatter agediff alder_atshock, absorb(cohort) line(none)
	graph export savings/moments_jeppe/atdeath/graphs/correlations/agediff-alderatshock.pdf, replace

	// Age difference and age at shock
	binscatter inheritance agediff, absorb(cohort) line(none)
	graph export savings/moments_jeppe/atdeath/graphs/correlations/agediff-alderatshock.pdf, replace
	
	
/* Lifecycle profiles */
local varstoplot networth_abs networth liqworth liqworth_abs illiqworth illiqworth_abs perminc inheritance inheritance_cens home_own
foreach var of local varstoplot {
	quietly {
		reg `var' b49.alder_atshock b2004.year b1960.cohort
		mat b = e(b)
		mat b = b[1,1..34]
		mat score _pred_`var' = b
		replace _pred_`var' = _pred_`var' + _b[_cons]
	}
}

collapse (count) count (mean) _pred* `varstoplot', by(alder_atshock cohort)

cap drop pickone
egen pickone = tag(alder_)

sort cohort alder

qui: levelsof cohort, local(levels)	
foreach var in `varstoplot' {
	local toplot
	foreach cohort of local levels {
		local toplot `toplot' (line `var' alder if cohort == `cohort', lc(black))
	}
	twoway `toplot' ///
					(line _pred_`var' alder if pickone, sort lw(medthick) lc(red)) ///
					, graphr(color(white)) legend(off) xtitle("Age") ytitle("") ///
					xlabel(25(5)60)
	graph export savings/moments_jeppe/atdeath/graphs/profiles/`var'.pdf, replace

}
preserve 
	keep alder _pred*
	duplicates drop
	sort alder
	rename _pred_* *
	export delimited using savings/toexport/profiles.csv, replace
restore

##################################
# Load estimates and export them #
##################################

foreach var in liqworth networth {
	clear all
	
	set obs 15
	gen time = _n + 94
	forval i = 95/109 {
		local real = `i'-100
		label define timef `i' "`real'", add
	}
	label val time timef

	estimates use savings/estimates/rev1/`var'_2_main
	
	gen beta_2 = .
	foreach i of numlist 95/98 99 100/109 {
		quietly {
			if `i'==99 {
				replace beta_2 = 0 if time == `i'
			}
			else{
				replace beta_2 = _b[`i'.timefid] if time == `i'
			}
		}
	}
	export delimited using savings/toexport/empirical_networth.csv, replace
}


