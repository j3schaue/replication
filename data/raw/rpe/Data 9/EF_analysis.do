capture log close /* close old log file */
set more off /* automatic scrolling */
clear
eststo clear
log using "EF.log", replace

/* START */

import delimited "EF_data.csv"

gen diff = wtamug - wtapen

gen log_wta_mug = ln(wtamug)
gen log_wta_pen = ln(wtapen)
gen logdiff = log_wta_mug - log_wta_pen

/* MAIN TEST */
ttest diff, by(hightreat) une
ttest logdiff, by(hightreat) une

/* MORE DATA*/

* average WTA
ttest wtamug, by(hightreat) une
ttest log_wta_mug, by(hightreat) une
ttest wtapen, by(hightreat) une
mean(wtamug)
mean(wtapen)

tab hightreat if wtamug==0
tab hightreat if wtapen==0

* Regressions 1-5 in Table 2

* Regression 1
eststo reg1 : regress log_wta_mug hightreat 

* Regression 2
* include a constant indicator for wtapen=9.57 (may be censored then)
gen pen_max_wta = 0
replace pen_max_wta = 1 if wtapen>9.5
eststo reg2 : regress log_wta_mug hightreat log_wta_pen pen_max_wta

* Regression 3 - add controls
eststo reg3 : regress log_wta_mug hightreat log_wta_pen female age pen_max_wta i.day

* Reg 4 is tobit that takes censoring at WTA_mug=9.57 into account
* ln(9.57) = 2.2586331
tobit log_wta_mug hightreat log_wta_pen female age pen_max_wta i.day, ul(2.2586331)
* display marginal effects
eststo reg4 : margins, dydx(*) atmeans post

* include cubic wtapen
eststo reg5 : regress log_wta_mug hightreat female age pen_max_wta c.wtapen##c.wtapen##c.wtapen i.day

esttab using EF_table2.rtf, replace title(Replication of Table II) label ///
	indicate("Day indicators = *.day" "Other controls = wtapen*") ///
	nomtitles noconstant se starlevels(* 0.10 ** 0.05 *** 0.01)
	
/* END */
log close
exit
