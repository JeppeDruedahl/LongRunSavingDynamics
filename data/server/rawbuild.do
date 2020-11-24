cd "D:\Data\workdata\705109"
clear all
set more off

*****************************
****** HOUSING DATASET ******
*****************************
use stataraw/grunddata/EJER1995, clear
gen year = 1995 -1
	forval year = 1996/2012 {
		append using stataraw/grunddata/EJER`year'
		replace year = `year' -1 if year==.
	}
egen bnr = group(ejdnr komnr)
bys pnr year: egen n_owned = count(bnr)
bys pnr year: egen pct_owned = total(ejerpct)

keep pnr year *_owned
duplicates drop

save data/housing, replace

******************************
***** ADDITIONAL  INCOME *****
******************************

use stataraw/grunddata/indk1995.dta, clear
gen year = 1995
forval year = 1996/2012 {
	append using stataraw/grunddata/indk`year'
	replace year = `year' if year==.
}
duplicates drop
save data/other_indk.dta, replace


*****************************
******* INCOME UPDATE *******
*****************************

use stataraw/grunddata/GRUND1995_NEW.dta, clear
gen year = 1995
forval year = 1996/2012 {
	append using stataraw/grunddata/GRUND`year'_NEW
	replace year = `year' if year==.
}
duplicates drop
rename erhv labor_income
rename ERHV labor_income_old
rename loenmv salary
egen pencontr = rowtotal(qarbpen qpripen)

/*
The next line of code drops a single line of data in the file, 
which represents an individual that we see twice in the data in 2009, 
one such line all filled with 0s except the personal ID. 
As this line of code technicaly constitutes microdata we cannot share it
and is therefore censored. 
The original dofile is kept in the server for replication purposes.
*/
* drop if XXXXXXX... //duplicate

save data/update1.dta, replace

******************************
******** MAIN DATASET ********
******************************
use stataraw/grunddata/GRUND1990, clear
gen year = 1990

forval year = 1991/2012 {
	append using stataraw/grunddata/GRUND`year'
	replace year = `year' if year==.
}

replace foed_dag = FOED_DAG if foed_dag == .
drop FOED_DAG

rename *, lower
duplicates drop
tab year

/*
The next line of code drops a single line of data in the file, 
which represents an individual that we see twice in the data in 2009, 
one such line all filled with 0s except the personal ID. 
As this line of code technicaly constitutes microdata we cannot share it
and is therefore censored. 
The original dofile is kept in the server for replication purposes.
*/
* drop if XXXXXXX... //duplicate


merge 1:1 pnr year using data/other_indk
drop if _merge==2
drop _merge


merge 1:1 pnr year using data/update1.dta
drop if _merge==2
drop _merge

/***** GENERATE VARIABLES *****/

/* XTSET */
	egen id = group(pnr)
	xtset id year

/* CPI normalization */
	gen cpi = .
	replace cpi =  66.09 if year == 1990
	replace cpi =  67.63 if year == 1991
	replace cpi =  69.12 if year == 1992
	replace cpi =  70.02 if year == 1993
	replace cpi =  71.41 if year == 1994
	replace cpi =  72.88 if year == 1995
	replace cpi =  74.43 if year == 1996
	replace cpi =  76.06 if year == 1997
	replace cpi =  77.45 if year == 1998
	replace cpi =  79.41 if year == 1999
	replace cpi =  81.70 if year == 2000
	replace cpi =  83.66 if year == 2001
	replace cpi =  85.62 if year == 2002
	replace cpi =  87.42 if year == 2003
	replace cpi =  88.48 if year == 2004
	replace cpi =  90.03 if year == 2005
	replace cpi =  91.75 if year == 2006
	replace cpi =  93.30 if year == 2007
	replace cpi =  96.49 if year == 2008
	replace cpi =  97.79 if year == 2009
	replace cpi = 100.00 if year == 2010
	replace cpi = 102.78 if year == 2011
	replace cpi = 105.23 if year == 2012

	local varstonorm brutto skatfriyd aktieindk skatmvialt_ny kurs* ///
		pant* obl* bank* koejd qpri* qarb* labor* salary resuink RESUINK_GL kappens ///
		
	foreach var of varlist `varstonorm' {
		replace `var' = `var'*100/cpi
	}
	
/* Income */

	egen dispinc = rowtotal(brutto skatfriyd aktieindk)
	replace dispinc = dispinc - skatmvialt_ny
	egen perminc = filter(dispinc), coef(45 25 15 10 5) lags(0/4) normalise

	keep if year>=1995
	
	*drop skatfriyd aktieindk skatmvialt_ny
	
/* Wealth */

	egen stocks = rowtotal(kursanp kursakt)
	replace stocks = max(0, stocks)
	egen bonds = rowtotal(pantakt oblakt)
	*gen bonds = 0
	replace bonds = max(0, bonds)
	rename bankakt liqassets
	replace liqassets = max(0, liqassets)
	rename bankgaeld liqdebts
	replace liqdebts = max(0, liqdebts)
	egen mortgage = rowtotal(oblgaeld pantgaeld)
	replace mortgage = max(0, mortgage)
	rename koejd hvalue
	replace hvalue = max(0, hvalue)
	gen hequity = hvalue - mortgage
	
	gen finw = stocks + bonds
	replace finw = max(0,finw)
	
	gen totassets = stocks + bonds + hvalue + liqassets
	gen totdebts = mortgage + liqdebts
	gen networth = totassets - totdebts
	
	drop kursanp kursakt stocks bonds pantgaeld
	

	/* Spouse */
	gen spo_pnr = efalle if civst == "G" | civst =="P"
	drop efalle
	
	/* Demographics */ 
	gen female = koen - 1
	
/***** DROP VARIABLES *****/
local varstodrop aegte_id civ_vfra van_vtil antboernf ///
	koen hf_* erhve* ansaar ansdage arbled jobkat persbr still tilk ///
	type arbfors arblh brutto haevpen haevpen korstoett korydi ///
	qaktiv* qpass* samlink skattot bel_ antpers doedsaars2 doedsaars3 qbisty
	

foreach var of varlist `varstodrop' {
	cap drop `var'
}

/***** HOUSING DATA *****/
merge 1:1 pnr year using data/housing
drop if _merge==2
drop _merge

replace n_owned = 0 if n_owned==. & year<2012
replace pct_owned = 0 if pct_owned==. & year<2012

/***** SAVE DATA *****/
compress
save data/core, replace


*********************************
********** SPOUSE DATA ********** 
*********************************
use spo_pnr year using data/core, clear

keep if spo_pnr!=""
duplicates drop

rename spo_pnr pnr
duplicates drop
merge 1:1 pnr year using data/core, keep(3) keepus(n_owned perminc dispinc pct_owned networth hequity liqassets finw liqdebts labor_income salary)
drop _merge
rename * spo_*
rename spo_year year

save data/spouse_vars, replace

**********************************
***** DEATHS AND INHERITANCE *****
**********************************
use data/core, clear
gen gross_bequest = max(0,max(0,hequity) + finw + liqassets - liqdebts)

cap drop sudden_death
gen sudden_death = inlist(substr(doedsaars1,1,3),"I21","I22","I23","I24","I46") ///
 | inlist(substr(doedsaars1,1,1), "V", "X", "Y") ///
 | inlist(substr(doedsaars1,1,3),"R96") if doedsaars1!=""

gen death_year = year(doddato)

keep pnr year gross_bequest sudden_death death_year
keep if inrange(death_year,1996,2012)

compress
save temp/grossbequests, replace

foreach p in m f {
	use stataraw/grunddata/FTDB2012, clear
	keep pnr pnr`p'
	rename pnr`p' pnr_sire
	save temp/fert_`p', replace
}
append using temp/fert_m
order pnr_sire pnr
duplicates drop
keep if pnr_sire!=""

bys pnr_sire: egen n_heirs = count(pnr!="")

expand 18
bys pnr_sire pnr: gen year = 1994+_n
tab year

rename pnr pnr_kid
rename pnr_sire pnr

merge m:1 pnr year using temp/grossbequests , keep(3)
drop _merge

bys pnr: egen lastyear = max(year)
gen inheritance = gross_bequest/n_heirs 
keep if year==lastyear
drop lastyear



/* for every kid, keep only largest inheritance*/
bys pnr_kid: egen large_inheritance = max(inheritance)
keep if large_inheritance == inheritance
drop year
drop large_inheritance

/* if parents have same wealth (98% 0), pick one at random (lower pnr) */
by pnr_kid: gen selected = pnr == pnr[1]
keep if selected 
drop selected

rename pnr sire_pnr
rename pnr_kid pnr

save data/inheritances, replace

**********************************
******* PARENTAL VARIABLES *******
**********************************
use sire_pnr using data/inheritances, clear
duplicates drop

rename sire_pnr pnr
merge 1:m pnr using data/core, keepus(year alder perminc n_owned pct_owned networth totassets hequity hvalue liqassets liqdebts)
keep if _merge==3
drop _merge

rename * sire_*
rename sire_year year
save data/sire_vars, replace

**********************************
******** MERGE FINAL DATA ********
**********************************

use data/core, clear
	/* spouse variables */
	merge m:1 spo_pnr year using data/spouse_vars
	drop _merge
	/* inheritance shock */
	merge m:1 pnr using data/inheritances
	keep if _merge==3
	drop _merge
	/* sire variables */
	merge m:1 sire_pnr year using data/sire_vars
	drop if _merge==2
	drop _merge

/* Remove taxes from inheritance */
	gen bundfr 		 = 184900*100/74.43 if death_year==1996
	replace bundfr = 186000*100/76.06 if death_year==1997
	replace bundfr = 191100*100/77.45 if death_year==1998
	replace bundfr = 196600*100/79.41 if death_year==1999
	replace bundfr = 203500*100/81.70 if death_year==2000
	replace bundfr = 210600*100/83.66 if death_year==2001
	replace bundfr = 216900*100/85.62 if death_year==2002
	replace bundfr = 224600*100/87.42 if death_year==2003
	replace bundfr = 231800*100/88.48 if death_year==2004
	replace bundfr = 236900*100/90.03 if death_year==2005
	replace bundfr = 242400*100/91.75 if death_year==2006
	replace bundfr = 248900*100/93.30 if death_year==2007
	replace bundfr = 255400*100/96.49 if death_year==2008
	replace bundfr = 264100*100/97.79 if death_year==2009
	replace bundfr = 264100 if death_year==2010
	replace bundfr = 264100*100/102.78 if death_year==2011
	replace bundfr = 264100*100/105.23 if death_year==2012
	
	replace inheritance = (inheritance*n_heirs-bundfr)*0.85 + bundfr if inheritance*n_heirs>bundfr
	replace inheritance = inheritance/n_heirs
	
/* some final variable management */
replace perminc = max(perminc, 60000)
replace spo_perminc = max(perminc, 60000)
foreach var in networth liqassets finw hequity hvalue liqdebts mortgage {
	replace `var' = `var'/perminc
}
gen x = perminc if (year+1)== death_year
bys id: egen perminc_base = total(x)
replace inheritance = inheritance/perminc_base
drop x

save data/workraw, replace


**********************************
******** SAMPLE SELECTION ********
**********************************
use data/workraw, clear
/* Working life */
	keep if inrange(alder,25,59)
	
	
sum inheritance if year==2000, d

/* Outliers in net worth */

	local count=0
	cap drop _sin*
	foreach var in networth {
		local ++count
		qui: sum `var', d
		gen _sin_`count' = !inrange(`var',r(p1),r(p99))
	}
	cap drop x
	egen x = rowtotal(_sin*)
	cap drop todrop
	bys id: egen todrop = max(x>0)
	drop x

	sum networth liqassets finw hequity hvalue liqdebts if !todrop

/* At least 3 years of observations before and after shock */
	gen timefromshock = year - death_year 
	bys pnr: egen minobs = min(timefromshock)
	bys pnr: egen maxobs = max(timefromshock)

/*Selection*/
	drop if todrop
	*keep if minobs<=-3 & maxobs>=3
	*drop minobs maxobs todrop

save data/workdata, replace




 
