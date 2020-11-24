clear all
set more off

/////////////////
// 1. load csv //
/////////////////

foreach model in $old_models {
	foreach suffix in "" $suffix_list {
		* a. standard
		import delimited model/csv/`model'`suffix'.txt, clear
		rename v1 age
		rename v5 inh
		rename v2 C_`model'`suffix'
		rename v3 A_`model'`suffix'
		rename v4 B_`model'`suffix'
		drop v6			// income - no need for now
		drop v7 		// permanent income - no need for now	
		destring age, replace force
		drop if age==.
		cap rename v8 type
		keep if age <= 59
		* generate id variable
		cap drop id
		gen id = floor((_n-1)/36)
		
		save data/local/datacustom/temp/`model'`suffix', replace
	}
}

foreach model in $new_models {
	foreach suffix in "" _unexp {
		* a. standard
		import delimited model/csv/`model'`suffix'.txt, clear
		rename v1 age
		rename v5 inh
		rename v2 C_`model'`suffix'
		rename v3 A_`model'`suffix'
		rename v4 B_`model'`suffix'
		drop v6			// income - no need for now
		drop v7 		// permanent income - no need for now	
		destring age, replace force
		drop if age==.
		cap rename v8 type
		keep if age <= 59
		* generate id variable
		cap drop id
		gen id = floor((_n-1)/36)
		
		save data/local/datacustom/temp/`model'`suffix', replace
	}
}


///////////////
// 2. append //
///////////////

local i = 1
foreach model in $old_models {
	foreach suffix in "" $suffix_list {
		if `i' == 1 {
			use using data/local/datacustom/temp/`model'`suffix', clear
			local i = 0
		}	
		else {
			merge 1:1 id age using data/local/datacustom/temp/`model'`suffix', keepus(C* A* B*)
			drop _merge
		}
	}
}
foreach model in $new_models {
	foreach suffix in "" _unexp {
		else {
			merge 1:1 id age using data/local/datacustom/temp/`model'`suffix', keepus(C* A* B*)
			drop _merge
		}
	}
}

//////////////////////
// 3. new variables //
//////////////////////

bys id (age): egen ageatdeath = max((inh!=0)*(age))
gen timefromshock = age - ageatdeath

foreach model in $old_models {
	foreach suffix in "" $suffix_list {
		gen NW_`model'`suffix' = A_`model'`suffix' + B_`model'`suffix'
	}
}

foreach model in $new_models {
	foreach suffix in "" _unexp {
		gen NW_`model'`suffix' = A_`model'`suffix' + B_`model'`suffix'
	}
}

/////////////
// 4. save //
/////////////

compress
if floor(c(version))==14 saveold data/local/datacustom/raw, version(13) replace
else save data/local/datacustom/raw, replace



***************************************
***** Generate lifecycle profiles *****
***************************************

/* networth */
// We want average, median, p25 and p75 of networth, for normal and nr 
use data/local/datacustom/raw, clear

*
keep age NW_* A_*

foreach model in $old_models {
	foreach suffix in "" $suffix_list {
		gen mean_NW_`model'`suffix' = NW_`model'`suffix'
		gen p50_NW_`model'`suffix' = NW_`model'`suffix'
		gen p25_NW_`model'`suffix' = NW_`model'`suffix'
		gen p75_NW_`model'`suffix' = NW_`model'`suffix'
	}
}
foreach model in $new_models {
	foreach suffix in "" _unexp {
		gen mean_A_`model'`suffix' = A_`model'`suffix'
		gen p50_A_`model'`suffix' = A_`model'`suffix'
		gen p25_A_`model'`suffix' = A_`model'`suffix'
		gen p75_A_`model'`suffix' = A_`model'`suffix'
	}
}

collapse (mean) mean_* (median) p50_* (p25) p25_* (p75) p75_* , by(age)


compress
if floor(c(version))==14 saveold data/local/datacustom/profiles, version(13) replace
else save data/local/datacustom/profiles, replace
