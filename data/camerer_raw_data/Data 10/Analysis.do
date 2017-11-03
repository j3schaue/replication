clear
set more off

use "Data.dta"

sum delegation if norec == 1 & type == 1
sum delegation if norec == 0 & type == 1

xi: probit delegation norec i.period if type==1, cluster(id)
mfx
drop _I*

table norec delegation if type == 1, content(mean earnings)

quietly {
		sum   earnings  if type == 1 & norec == 0 & delegation == 0
		local a0_d0 = r(mean)
		sum   earnings  if type == 1 & norec == 0 & delegation == 1
		local a0_d1 = r(mean)
		sum   earnings  if type == 1 & norec == 1 & delegation == 0
		local a1_d0 = r(mean)
		sum   earnings  if type == 1 & norec == 1 & delegation == 1
		local a1_d1 = r(mean)
}

* Actual Return in PHigh25
dis `a0_d1'/`a0_d0'
* Actual Return in High NoRec
dis `a1_d1'/`a1_d0'
