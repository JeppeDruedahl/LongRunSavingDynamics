
cap program drop timerec
program define timerec
version 11

syntax , SAVEfile(string) Timer(integer) [TIMENote(string)]
	
	
**time recording
timer list
local exph1 = floor( r(t`timer')/3600  )
local expm1 = floor( r(t`timer')/60  ) - `exph1'*60
local exps1 = floor( r(t`timer')  ) - `expm1'*60 - `exph1'*3600

di "Run time: `exph1' hour(s), `expm1' min(s) and `exps1' sec(s)"

local today = c(current_date)
*if length("`today'")==11 local today = "`today'"
*else local today = " `today'"

local skd=11-length("`today'")
local skh=6-length("`exph1'")
local skm=2-length("`expm1'") 
local sks=2-length("`exps1'")

cap mkdir timerec
cap confirm file `savefile'
local _err = _rc!=0
cap file close mylog
if `_err'==0 {
	file open mylog using `savefile',read
		local _linenum = 1
		file read mylog line1
		while r(eof)==0 {
			local _linenum= `_linenum'+1
			file read mylog line`_linenum'
		}
	file close mylog
}
file open mylog using `savefile', write replace
file write mylog _skip(`skd')"`today':" _skip(`skh') "`exph1' hour(s), " _skip(`skm')"`expm1' min(s) and " _skip(`sks')"`exps1' sec(s)" _skip(3)"`timenote'" _n
if `_err'==0 {
	forval i=1(1)`_linenum' {
		file write mylog "`line`i''" _n
	}
}
file close mylog

end
