clear all
set more off
use datacustom/raw, clear
destring *, replace
keep if age>=25
drop *_nr *_nvr

//////////////////
// 1. adj. time //
//////////////////

replace timefromshock = timefromshock + 100
replace timefrom = timefrom if timefrom>=100

cap drop timefid
gen timefid = (timefrom>=96)*timefrom

gen timefid2 = timefrom
replace timefid2 = 99 if timefid2==80

////////////////////
// 2. regressions //
////////////////////

keep if age < 55
keep if inrange(timefid, 80, 109)
bys timefid age: egen aw = count(timefid)
replace aw = _N/aw

/* Full sample*/ 
foreach model in $old_models $new_models {
	foreach suffix in "" "_unexp" {
		foreach var in A B NW {
			
			reg `var'_`model'`suffix' i.age i.age b99.timefid [aw=aw]
			est save data/local/estimates/simulated/`var'_`model'`suffix'_all, replace
			
		}
	}
}
*


