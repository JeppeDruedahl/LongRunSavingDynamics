cap program drop mkpath
program define mkpath
version 13

syntax anything

foreach p in `anything' {
	gettoken first p : p, parse("/\")
	cap mkdir "`first'"
	
	tokenize `p', parse("/\")
	local pl
	local i=1
	while "``i''"!="" {
		if !inlist("``i''","/","\") {
			local pl "`pl'/``i''"
			cap mkdir "`first'/`pl'"
		}
		local i=`i'+1
	}
}
end
