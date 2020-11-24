**********************************************
** GENERATE FULL DEATH SAMPLE FOR COMPUTING **
** MARGINAL PROBABILITIES OF DEATH BY AGE		**
**********************************************
cd "D:\Data\workdata\705109"
clear all
set more off



forval y = 1980/2010 {
	use pnr alder /*qaktivf qpassiv*/ doddato doedsaars1ALESSANDRO using stataraw/grunddata/GRUND`y'
	gen year = `y'
	save temp/deathstuff/death_`y', replace
}
use temp/deathstuff/death_1980
forval y = 1981/2010 {
	append using temp/deathstuff/death_`y'
}
save temp/rawdeath, replace
*
keep if doddato!=.
gen cohort = year - alder


rename doedsaars1ALESSANDRO doedsaars1ALESSANDRO
cap drop sudden_death
gen sudden_death = inlist(substr(doedsaars1ALESSANDRO,1,3),"I21","I22","I23","I24","I46") ///
 | inlist(substr(doedsaars1ALESSANDRO,1,1), "V", "X", "Y") ///
 | inlist(substr(doedsaars1ALESSANDRO,1,3),"R96") if doedsaars1ALESSANDRO!=""

gen death_year = year(doddato)

gen alderatdeath = death_year - cohort

keep pnr alderatdeath death_year sudden cohort 
duplicates drop

save temp/workdeath, replace








