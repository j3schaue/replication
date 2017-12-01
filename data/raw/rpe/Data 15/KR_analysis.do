cap log close
set memory 500m
clear
log using "main.log", replace
/* START */

/* Load full dataset */
use "KR_data.dta", clear

drop if treatment == 1 /* drop instructions */
drop if treatment == 17 /* drop mid game instructions */
drop if treatment == 33 /* drop last round announcement */
drop if treatment == 35 /* drop payment round */
by treatment, sort: gen play = _n == 1 /* "play" is the round of the game */
replace play = sum(play)

* Donation decisions are only made once per round so we need to drop all
* but the first period even when looking at averages since otherwise the length
* of the round will affect the average.
drop if period ~= 1 /* we are just running regression on first period of each round */

*bysort lifesharers: sum donate if play == 1
*bysort lifesharers: sum donate if play == 16

table play lifesharers, c(m donate) /* TABLE 1*/

keep donate lifesharers u_subject period play

* Test comparison to regression 1 (all 1-31 rounds)
reg donate lifesharers, robust cluster(u_subject)
est tab, p(%12.10g)

* Run LPM (no controls needed) <-- MAIN REPLICATION TEST
drop if play >= 16 /* and just on 15 first games */
reg donate lifesharers, robust cluster(u_subject)

/* COMPARE TO ORIGINAL KESSLER & ROTH RE-ESTIMATED */
clear
use "KR_original_data.dta"
drop if period ~= 1
table play game, c(m donate)
* drop if rebate==1 | decreasecost==1 | rebatefirst==1 | decreasefirst==1 | rebatesecond == 1 | decreasesecond == 1
* we should keep all observations that are not in treatment to get same data as original paper
drop if rebate==1 | decreasecost==1
* regression 1 equivalent (should be 0.306 effect size)
xi: reg donate lifesharers, robust cluster(u_subject)
est tab, p(%12.10g)
* only first periods (gives new estimate 0.383)
drop if play>=16
xi: reg donate lifesharers, robust cluster(u_subject)

/* END */
log close
exit
