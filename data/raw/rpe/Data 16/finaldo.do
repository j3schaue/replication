cd "x"


use t1_m2, clear
keep if Period==0
local i=1
while `i' <= 2{
	local j=2
	while `j' <= 12{
	capture append using t`i'_m`j'
	local j = `j' + 1
}
local i = `i' + 1
}
save "contracts_all", replace



******************************************
//Step2) extract trades from contracts table
******************************************

use contracts_all, clear
keep if Seller >0 & Buyer >0

**Create a table containing FV-information and merge it.
sort treatment Period
merge treatment Period using fv
drop _merge
save trades_all, replace




**************************************
//Calcualtion of BUBBLE MEASURES
**************************************
use trades_all, clear
egen volume=sum(q), by(treatment markt Period)
gen p2=p*q
egen sum_p=sum(p2), by(treatment markt Period)
gen ampreis=sum_p/volume

collapse ampreis fv volume, by(treatment markt Period)
gen diff=ampreis-fv
gen abs_diff=abs(diff)

egen mean_fv=mean(fv), by(treatment markt)

**RAD**
gen rad_help=abs_diff/mean_fv
egen rad=mean(rad_help), by(treatment markt)

**RD**
gen rd_help=diff/mean_fv
egen rd=mean(rd_help), by(treatment markt)

collapse rad rd, by(treatment markt)
sort treatment markt

save measures_data, replace



******************************
//Run Mann-Whitney-tests
******************************

use measures_data, clear

ttest rad, by(treatment)
ranksum rad, by(treatment)

ttest rd, by(treatment)
ranksum rd, by(treatment)

