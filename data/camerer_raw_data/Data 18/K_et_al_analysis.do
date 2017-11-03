capture log close /* close old log file */
set more off /* automatic scrolling */
clear
eststo clear
log using "K_et_al.log", replace
set seed 1 //Set seed for comparability

/* START */
use K_et_al_data.dta, clear

* drop practice round
drop if round==1

*The dependent variable is 1 if the participant has chosen to give money to a player with lower rank and 0 otherwise.
generate givetolower = 0
replace givetolower=1 if decision=="give to lower rank"

* Give to lower should depend on choice to give to lower in choice set, not lower relative rank.
* So rank 1 giving to 2 is not givetolower.
replace givetolower=0 if rank==1 & choice==2
* However when the person placed last gives to the 5th ranked it is givetolower.
replace givetolower=1 if rank==6 & choice==5

*Main indep variable
gen secondfromlast = 0
replace secondfromlast=1 if rank==5

* Dummies for first/last place
gen lastplace = 0
replace lastplace=1 if rank==6
gen firstplace = 0
replace firstplace=1 if rank==1

* There should be no missing data but test just in case
replace givetolower=. if decision=="NULL"
replace givetolower=. if missing(decision)
replace secondfromlast=. if missing(rank)
replace lastplace=. if missing(rank)
replace firstplace=. if missing(rank)

keep player_id game_id round rank choice givetolower secondfromlast firstplace lastplace

* Also create second or third from last for robustness checks
generate secondorthird_fromlast = secondfromlast
replace secondorthird_fromlast = 1 if rank==4

/* RUN ANALYSIS */

* GRAPH (Figure 4)
* "ib5" to use rank 5 as base indicator
probit givetolower ib5.rank, cluster(player_id)
margins rank, atmeans

* REGRESSIONS

* for regression 1 in table II
eststo clear
probit givetolower secondfromlast firstplace lastplace i.game_id i.round, cluster(player_id)
eststo reg1: margins, atmeans dydx(*) post

* same with LPM to avoid incidental-parameters problem
* eststo reg2: qui regress givetolower secondfromlast firstplace lastplace i.game_id i.round, robust cluster(player_id)

* regression 2 in table II
probit givetolower secondorthird_fromlast firstplace lastplace i.game_id i.round, cluster(player_id)
eststo reg3: margins, atmeans dydx(*) post

esttab using K_et_altable2.rtf, replace title(Replication of Table II) label ///
	indicate("Game f.e. = *.game_id" "Round f.e. = *.round") ///
	mtitles("Probit 1" "LPM" "Probit 3") ///
	drop(firstplace lastplace) ///
	compress noconstant se starlevels(* 0.10 ** 0.05 *** 0.01)

esttab, replace title(Replication of Table II) label ///
	indicate("Game f.e. = *.game_id" "Round f.e. = *.round") ///
	mtitles("Probit 1" "LPM" "Probit 3") ///
	drop(firstplace lastplace) ///
	compress plain noconstant p starlevels(* 0.10 ** 0.05 *** 0.01)


/* ROBUSTNESS */
/* TABLE 7 (APPENDIX)*/

eststo clear

* Robustness check with questionaire data (recreation of table 7 in appendix)

* Original version but with 2nd or 3rd place as explanatory var.
qui probit givetolower secondorthird_fromlast firstplace lastplace i.game_id i.round, cluster(player_id)
eststo reg1: qui margins, atmeans dydx(*) post

merge m:1 player_id using "K_et_al_questionaire.dta", keep(master match)
drop if _merge==1
drop if missing(male)
drop if missing(black)
drop if missing(hispanic)

*Check change to original spec. from dropping data with no questionaire answers
qui probit givetolower secondorthird_fromlast firstplace lastplace i.game_id i.round, cluster(player_id)
eststo reg2: qui margins, atmeans dydx(*) post

*Add demographic controls
qui probit givetolower secondorthird_fromlast male black hispanic age political_views religiosity firstplace lastplace i.game_id i.round, cluster(player_id)
eststo reg3: qui margins, atmeans dydx(*) post

*with political views interaction
*qui probit givetolower i.secondorthird_fromlast##c.political_views male black hispanic age religiosity firstplace lastplace i.game_id i.round, cluster(player_id)
*eststo reg4: qui margins, atmeans dydx(*) post

*with religiosity interaction
*qui probit givetolower i.secondorthird_fromlast##c.religiosity male black hispanic age political_views firstplace lastplace i.game_id i.round, cluster(player_id)
*eststo reg5: qui margins, atmeans dydx(*) post

*LPM
eststo reg6: qui regress givetolower secondorthird_fromlast male black hispanic age political_views religiosity firstplace lastplace i.game_id i.round i.player_id, cluster(player_id)

esttab using K_et_altable7.rtf, replace title(Replication of Appendix Table 7) label ///
	indicate("Game f.e. = *.game_id" "Round f.e. = *.round" "Indiv. f.e. = *.player_id") ///
	mtitles("Probit 1" "Probit 2" "Probit 3" "LPM") ///
	compress noconstant se starlevels(* 0.10 ** 0.05 *** 0.01)
	
esttab, replace title(Replication of Appendix Table 7) label ///
	indicate("Game f.e. = *.game_id" "Round f.e. = *.round" "Indiv. f.e. = *.player_id") ///
	mtitles("Probit 1" "Probit 2" "Probit 3" "LPM") ///
	compress plain noconstant p starlevels(* 0.10 ** 0.05 *** 0.01) 
*/

/* END */
log close
exit
