clear
set more off

use "Data.dta"

* Dummy for Observations gathered in September/October 2015 (seoc = 1)
* If seoc = 0, Observations were gathered in June/July 2015
gen     seoc = 0
replace seoc = 1 if sid > 96

sum     profit_main if treat_hi == 1, det
sum     profit_main if treat_hi == 0, det

ttest   temp, by(seoc)

* Hypothesis Test for Replication Result
reg     profit_main treat_hi
reg     profit_main treat_hi      if seoc == 0
reg     profit_main treat_hi      if seoc == 1

* Controlling for Room Temerature
reg     profit_main treat_hi temp
reg     profit_main treat_hi temp if seoc == 0
reg     profit_main treat_hi temp if seoc == 1

* [Temperature x Treatment]
gen txt = treat_hi * temp
reg     profit_main treat_hi temp txt
* [Sept/Oct x Treatment]
gen sxt = seoc * treat_hi
reg     profit_main treat_hi seoc sxt
