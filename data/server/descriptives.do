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



*****************************
******** SET GLOBALS ********
*****************************
global outcomes networth liqworth illiqworth liqassets liqdebts finw hequity hvalue mortgage
global outcomes_abs	networth_abs liqworth_abs illiqworth_abs liqassets_abs liqdebts_abs finw_abs hequity_abs hvalue_abs mortgage_abs
global other_outcomes n_owned home_own multi_own married antboernh spo_networth spo_n_owned dispinc labor_income salary has_labor has_qpripen qpripen_abs qpripen_norm has_qarbpen qarbpen_abs qarbpen_norm
global household_outcomes hh_networth hh_dispinc hh_labor_income
global mainfes i.year#i.cohort


******************************
*** CREATE ADDITIONAL VARS ***
******************************

/* Asset grouping */ 
gen liqworth = liqassets - liqdebts
gen illiqworth = hequity + finw 


/* Temporary - get back to base values*/
foreach var in $outcomes inheritance {
	gen `var'_abs = `var'*perminc/1000
}
replace dispinc = dispinc/1000
replace labor_income = labor_income/1000
replace salary = salary/1000

gen 		inhgroup = inheritance>1/12
replace inhgroup = 2 if inheritance>1


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


gen positive_b = illiqworth>0

cap drop wealth_before
cap drop countobs

/*
bys id: egen wealth_before = total(inrange(timefromshock,97,99)*networth)
by id: egen countobs = total(inrange(timefromshock,97,99))
replace wealth_before = wealth_before/countobs
sum wealth_before if timefromshock==99, d
cap drop high_wealth
gen high_wealth = wealth_before>r(p50)
cap drop wealth_q
gen wealth_q = 1 if wealth_before<=r(p25)
replace wealth_q = 2 if wealth_before>r(p25) & wealth_before<=r(p50)
replace wealth_q = 3 if wealth_before>r(p50) & wealth_before<=r(p75)
replace wealth_q = 4 if wealth_before>r(p75)
*/

bys id (year): egen _x = mean(networth) if inrange(timefid,96,99)
forval i = 0/1 {
	cap drop high_wealth
	cap drop _x
	bys alder: egen _x = median(networth) if timefid==99
	bys id (year): egen high_wealth = max(networth>_x)
	cap drop _x
}

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


replace spo_networth = spo_networth/spo_perminc
replace spo_hequity = spo_hequity/spo_perminc

replace spo_perminc = 0 if spo_pnr ==""

foreach var in networth {
	replace spo_`var' = 0 if spo_pnr ==""
	cap drop hh_`var'
	gen hh_`var' = (`var'*perminc + spo_`var'*spo_perminc)/(perminc + spo_perminc)
}
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

replace perminc = perminc/1000

cap drop inheritance_cens
gen inheritance_cens = min(inheritance, 15) if inheritance!=.

/**********************************/
****** HISTOGRAMS INHERITANCE ******
/**********************************/
cap drop _x _y
kdensity inheritance_cens if sudden==1, n(75) gen(_x _y) bw(0.15) nograph

twoway 	(line _y _x, color(black) sort) ///
	if _x>1, graphr(color(white)) ///
	xtitle("Inheritance (over permanent income, censored at 15)") ytitle("Density") ///
	xline(1, lp(-) lc(black))
graph export savings\graphs\rf_levels/inheritance_distribution.pdf, replace


/********************************/
******* DESCRIPTIVE TABLE *******
/********************************/
mkpath savings\tabs/${cdate}
keep if timefromshock==99
/*
GROUPS
sudden:			1 Unexpected		0 Expected
inhrgroup:	0 Placebo				2 Treatment
old: 				0 <= 40 yo			1 > 40 yo
*/

local varstodescribe inheritance_abs perminc inheritance networth liqassets ///
		liqdebts finw hequity hvalue mortgage home_own multi_own dispinc married ///
		death_year alder_atshock sire_alder_atshock 
local n_rows = wordcount("`varstodescribe'")

mat means = J(`n_rows',3,.)
mat append = J(1,3,.)
local count = 0
foreach var of local varstodescribe {
	local ++count
	qui: sum `var', meanonly
	mat means[`count',1] = r(mean)
	qui: sum `var', meanonly, if sudden==1
	mat means[`count',2] = r(mean)
	qui: sum `var', meanonly, if sudden==1 & inhgroup==2
	mat means[`count',3] = r(mean)
	
	/*qui: sum `var', meanonly, if sudden==1 & inhgroup==2 & old==0
	mat means[`count',4] = r(mean)
	qui: sum `var', meanonly, if sudden==1 & inhgroup==2 & old==1
	mat means[`count',5] = r(mean)
	*/
}
qui: count
mat append[1,1] = r(N)
qui: count if sudden==1
mat append[1,2] = r(N)
qui: count if sudden==1 & inhgroup==2
mat append[1,3] = r(N)
/*
qui: count if sudden==1 & inhgroup==2 &  old==0
mat append[1,4] = r(N)
qui: count if sudden==1 & inhgroup==2 &  old==1
mat append[1,5] = r(N)
*/

local varlabels "Imputed inheritance, 1000DKK" "Permanent income, 1000DKK" "Imputed inheritance, normalized" "Net worth" ///
	"$\quad-$ Liquid assets"  "$\quad-$ Uncollateralized debts" "$\quad-$ Financial investments" "$\quad-$ Housing equity" ///
	"$\qquad-$ Housing value" "$\qquad-$ Mortgage" "$\qquad-$ Home owner" "$\qquad-$ Owner of 2+ units" ///
	"Disposable income" "Married" "Year of inheritance" "Age at inheritance" "Parental age at death" 	
write_mats means using savings\tabs/descriptives.txt, ///
		names(`""`varlabels'""') space(0.5) append(append) appendnames(`""\# episodes""')


/*
Distribution of wealth and income across groups
*/
keep if sudden==1

twoway 	(kdensity networth if inhgroup==2, lc(black)) ///
		(kdensity networth if inhgroup!=2, lc(red) lp(-)) ///
		, graphr(color(white)) ///
		xtitle("Net worth, normalized by permanent income") ytitle("") ///
		legend(region(color(white)) order(1 "Large potential inheritance" 2 "Remaining sample"))
graph export savings\graphs\descriptives\distribution_networth.pdf, replace 
		

twoway 	(kdensity perminc if inhgroup==2, lc(black)) ///
		(kdensity perminc if inhgroup!=2, lc(red) lp(-)) ///
		, graphr(color(white)) ///
		xtitle("Permanent income (1000 DKK)") ytitle("") ///
		legend(region(color(white)) order(1 "Large potential inheritance" 2 "Remaining sample"))
graph export savings\graphs\descriptives\distribution_perminc.pdf, replace 


