program define cdate
	local date = c(current_date)
	local date = date("`date'","DMY")
	local y =substr(string(year(`date')),3,2)
	local m = cond(length(string(month(`date')))==1,string(0)+string(month(`date')),string(month(`date')))
	local d = cond(length(string(day(`date')))==1,string(0)+string(day(`date')),string(day(`date')))
	global cdate = "`y'"+"`m'"+"`d'"
end
