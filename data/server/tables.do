cd "D:\Data\workdata\705109"
clear all
set more off
adopath++ "ado\"
cdate

**************************************
***** REDUCED  FORM - CATEGORIES *****
**************************************
/*
One table for each wealth component (appendix)
Show only FE?
*/


mkpath savings\tabs/${cdate}



/***** MAIN TABLE ******/
/* Effect before 1 (99.) and after 1 (101.), 5 (105.), 9(109.) years */
	/* Wealth - main */
	local varstotab networth liqassets hequity finw liqdebts  
	local n_rows = wordcount("`varstotab'")
	mat B = J(`n_rows',8,.)
	mat SE = J(`n_rows',8,.)
	local countrow = 0
	foreach var of local varstotab {
		local ++countrow
		estimates use savings/estimates/rev1/`var'_abs_2_main
		
			mat B[`countrow',1] = _b[98.timefid]
			mat SE[`countrow',1] = _se[98.timefid]
			mat B[`countrow',2] = _b[101.timefid]
			mat SE[`countrow',2] = _se[101.timefid]
			mat B[`countrow',3] = _b[105.timefid]
			mat SE[`countrow',3] = _se[105.timefid]
			mat B[`countrow',4] = _b[109.timefid]
			mat SE[`countrow',4] = _se[109.timefid]
			
			estimates use savings/estimates/rev1/`var'_2_main
			mat B[`countrow',5] = _b[98.timefid]
			mat SE[`countrow',5] = _se[98.timefid]
			mat B[`countrow',6] = _b[101.timefid]
			mat SE[`countrow',6] = _se[101.timefid]
			mat B[`countrow',7] = _b[105.timefid]
			mat SE[`countrow',7] = _se[105.timefid]
			mat B[`countrow',8] = _b[109.timefid]
			mat SE[`countrow',8] = _se[109.timefid]
			
	}
	* print
		write_mats B SE using savings\tabs/empirics_2.txt, ///
			 space(0.5) /*stars(B SE)*/ ///
			names(`""Net worth" "$\quad-$ Liq. assets" "$\quad-$ Housing equity" "$\quad-$ Fin. investments" "$\quad-$ Unc. debts""')
	*
	/* Wealth - placebo */

	local varstotab networth liqassets hequity finw liqdebts  
	local n_rows = wordcount("`varstotab'")
	mat B = J(`n_rows',8,.)
	mat SE = J(`n_rows',8,.)
	local countrow = 0
	foreach var of local varstotab {
		local ++countrow
		estimates use savings/estimates/rev1/`var'_abs_0_main
		
			mat B[`countrow',1] = _b[98.timefid]
			mat SE[`countrow',1] = _se[98.timefid]
			mat B[`countrow',2] = _b[101.timefid]
			mat SE[`countrow',2] = _se[101.timefid]
			mat B[`countrow',3] = _b[105.timefid]
			mat SE[`countrow',3] = _se[105.timefid]
			mat B[`countrow',4] = _b[109.timefid]
			mat SE[`countrow',4] = _se[109.timefid]
			
			estimates use savings/estimates/rev1/`var'_0_main
			mat B[`countrow',5] = _b[98.timefid]
			mat SE[`countrow',5] = _se[98.timefid]
			mat B[`countrow',6] = _b[101.timefid]
			mat SE[`countrow',6] = _se[101.timefid]
			mat B[`countrow',7] = _b[105.timefid]
			mat SE[`countrow',7] = _se[105.timefid]
			mat B[`countrow',8] = _b[109.timefid]
			mat SE[`countrow',8] = _se[109.timefid]
			
	}
	* print
		write_mats B SE using savings\tabs/empirics_0.txt, ///
			 space(0.5) /*stars(B SE)*/ ///
			names(`""Net worth" "$\quad-$ Liq. assets" "$\quad-$ Housing equity" "$\quad-$ Fin. investments" "$\quad-$ Unc. debts""') 
	*	*
	
	/* OTHER VARS */
	/* income and labor supply */
	local varstotab dispinc has_labor labor_income salary spo_dispinc spo_has_labor spo_labor_income spo_salary
	
	local n_rows = wordcount("`varstotab'")
	mat B = J(`n_rows',4,.)
	mat SE = J(`n_rows',4,.)
	local countrow = 0
	foreach var of local varstotab {
		local ++countrow
			estimates use savings/estimates/rev1/`var'_2_main
			
			mat B[`countrow',1] = _b[98.timefid]
			mat SE[`countrow',1] = _se[98.timefid]
			mat B[`countrow',2] = _b[101.timefid]
			mat SE[`countrow',2] = _se[101.timefid]
			mat B[`countrow',3] = _b[105.timefid]
			mat SE[`countrow',3] = _se[105.timefid]
			mat B[`countrow',4] = _b[109.timefid]
			mat SE[`countrow',4] = _se[109.timefid]
	}
	* print
	write_mats B SE using savings\tabs/other_2_inc.txt, ///
		 space(0.5) /*stars(B SE)*/ ///
		names(`""Net disposable income" "Has earnings" "Gross earnings" "Gross salary" "Spouse net disposable income" "Spouse has earnings" "Spouse gross earnings" "Spouse gross salary" "')

	/* pension contributions */
	local varstotab qarbpen_norm qpripen_norm
	
	local n_rows = wordcount("`varstotab'")
	mat B = J(`n_rows',4,.)
	mat SE = J(`n_rows',4,.)
	local countrow = 0
	foreach var of local varstotab {
		local ++countrow
			estimates use savings/estimates/rev1/`var'_2_main
			
			mat B[`countrow',1] = _b[98.timefid]
			mat SE[`countrow',1] = _se[98.timefid]
			mat B[`countrow',2] = _b[101.timefid]
			mat SE[`countrow',2] = _se[101.timefid]
			mat B[`countrow',3] = _b[105.timefid]
			mat SE[`countrow',3] = _se[105.timefid]
			mat B[`countrow',4] = _b[109.timefid]
			mat SE[`countrow',4] = _se[109.timefid]
	}
	* print
		write_mats B SE using savings\tabs/other_2_pens.txt, ///
			 space(0.5) /*stars(B SE)*/ ///
			names(`""Employment scheme" "Personal funds""')
			
	/* household*/
	local varstotab married antboernh spo_networth hh_networth
	
	local n_rows = wordcount("`varstotab'")
	mat B = J(`n_rows',4,.)
	mat SE = J(`n_rows',4,.)
	local countrow = 0
	foreach var of local varstotab {
		local ++countrow
			estimates use savings/estimates/rev1/`var'_2_main
			
			mat B[`countrow',1] = _b[98.timefid]
			mat SE[`countrow',1] = _se[98.timefid]
			mat B[`countrow',2] = _b[101.timefid]
			mat SE[`countrow',2] = _se[101.timefid]
			mat B[`countrow',3] = _b[105.timefid]
			mat SE[`countrow',3] = _se[105.timefid]
			mat B[`countrow',4] = _b[109.timefid]
			mat SE[`countrow',4] = _se[109.timefid]
	}
	* print
		write_mats B SE using savings\tabs/other_2_hhd.txt, ///
			 space(0.5) /*stars(B SE)*/ ///
			names(`""Married" "\# children" "Spouse net worth$^{a}$" "Household net worth$^{b}$" "') 

			
	/* HOUSING - treated */
	local varstotab hequity hvalue home_own multi_own mortgage 
	local n_rows = wordcount("`varstotab'")
	mat B = J(`n_rows',4,.)
	mat SE = J(`n_rows',4,.)
	local countrow = 0
	foreach var of local varstotab {
		local ++countrow
			estimates use savings/estimates/rev1/`var'_2_main
			
			mat B[`countrow',1] = _b[98.timefid]
			mat SE[`countrow',1] = _se[98.timefid]
			mat B[`countrow',2] = _b[101.timefid]
			mat SE[`countrow',2] = _se[101.timefid]
			mat B[`countrow',3] = _b[105.timefid]
			mat SE[`countrow',3] = _se[105.timefid]
			mat B[`countrow',4] = _b[109.timefid]
			mat SE[`countrow',4] = _se[109.timefid]
	}
	* print
		write_mats B SE using savings\tabs/housing_2.txt, ///
			 space(0.5) /*stars(B SE)*/ ///
			names(`""Housing equity" "$\quad-$ Housing value" "$\qquad-$ Home owner" "$\qquad-$ Multiple real estate" "$\quad-$ Mortgage" "') 
	*	
	
	/* ROLE OF LQ COSTRAINTS  */
local varstotab networth liqassets hequity finw liqdebts  
	local n_rows = wordcount("`varstotab'")
	mat B = J(`n_rows',8,.)
	mat SE = J(`n_rows',8,.)
	local countrow = 0
	foreach var of local varstotab {
		local ++countrow
		estimates use savings/estimates/rev1/`var'_2_uncostr
		
			mat B[`countrow',1] = _b[98.timefid]
			mat SE[`countrow',1] = _se[98.timefid]
			mat B[`countrow',2] = _b[101.timefid]
			mat SE[`countrow',2] = _se[101.timefid]
			mat B[`countrow',3] = _b[105.timefid]
			mat SE[`countrow',3] = _se[105.timefid]
			mat B[`countrow',4] = _b[109.timefid]
			mat SE[`countrow',4] = _se[109.timefid]
			
			estimates use savings/estimates/rev1/`var'_2_costr
			mat B[`countrow',5] = _b[98.timefid]
			mat SE[`countrow',5] = _se[98.timefid]
			mat B[`countrow',6] = _b[101.timefid]
			mat SE[`countrow',6] = _se[101.timefid]
			mat B[`countrow',7] = _b[105.timefid]
			mat SE[`countrow',7] = _se[105.timefid]
			mat B[`countrow',8] = _b[109.timefid]
			mat SE[`countrow',8] = _se[109.timefid]
			
	}
	* print
		write_mats B SE using savings\tabs/role_lqconstr_2.txt, ///
			 space(0.5) /*stars(B SE)*/ ///
			names(`""Net worth" "$\quad-$ Liq. assets" "$\quad-$ Housing equity" "$\quad-$ Fin. investments" "$\quad-$ Unc. debts""')

/***** APPENDIX - networth and liqworh, FN and event *****/
local esttotabs networth_event_yearcohort networth_FN_yearcohort_nobalance networth_FN_yearcohort_balance liqassets_event_yearcohort liqassets_FN_yearcohort_nobalance liqassets_FN_yearcohort_balance

local n_cols = wordcount("`esttotabs'")
	mat B = J(11,`n_cols',.)
	mat SE = J(11,`n_cols',.)
	local countcol = 0
	foreach est of local esttotabs {
		local ++countcol
		est use savings/estimates/FN/`est'
		
		local countrow = 0
		foreach t of numlist 97/98 100/108 {
			local ++countrow 
			if strmatch("`est'","*_FN_*") == 0 {
				mat B[`countrow',`countcol'] = _b[`t'.timefid]
				mat SE[`countrow',`countcol'] = _se[`t'.timefid]
			}
			else if strmatch("`est'","*_FN_*") == 1 {
				mat B[`countrow',`countcol'] = _b[1.treat#`t'.timefid]
				mat SE[`countrow',`countcol'] = _se[1.treat#`t'.timefid]
			}
		}
	}
	write_mats B SE using savings\tabs/FN_VS_event.txt, ///
			 space(0.5) stars(B SE) ///
			names(`""-3" "-2" "0" "1" "2" "3" "4" "5" "6" "7" "8""')

/***** APPENDIX - all estimated coeff *****/
local main networth_2_main liqassets_2_main hequity_2_main finw_2_main liqdebts_2_main 
local main_abs networth_abs_2_main liqassets_abs_2_main hequity_abs_2_main finw_abs_2_main liqdebts_abs_2_main 
local placebo networth_0_main liqassets_0_main hequity_0_main finw_0_main liqdebts_0_main 
local placebo_abs networth_abs_0_main liqassets_abs_0_main hequity_abs_0_main finw_abs_0_main liqdebts_abs_0_main 
local housing hequity_2_main hvalue_2_main home_own_2_main multi_own_2_main mortgage_2_main
local income_pension dispinc_2_main labor_income_2_main salary_2_main qarbpen_norm_2_main qpripen_norm_2_main 
local household antboernh_2_main married_2_main spo_networth_2_main hh_networth_2_main

foreach list in main main_abs placebo placebo_abs housing income_pension houshold {
	local esttotabs ``list''
	
	local n_cols = wordcount("`esttotabs'")
	mat B = J(14,`n_cols',.)
	mat SE = J(14,`n_cols',.)
	local countcol = 0
	foreach est of local esttotabs {
		local ++countcol
		est use savings/estimates/rev1/`est'
		
		local countrow = 0
		foreach t of numlist 95/98 100/109 {
			local ++countrow 
			mat B[`countrow',`countcol'] = _b[`t'.timefid]
			mat SE[`countrow',`countcol'] = _se[`t'.timefid]
		}
	}
	write_mats B SE using savings\tabs/`list'.txt, ///
			 space(0.5) stars(B SE) ///
			names(`""-5" "-4" "-3" "-2" "0" "1" "2" "3" "4" "5" "6" "7" "8" "9""')

}	
	
	
