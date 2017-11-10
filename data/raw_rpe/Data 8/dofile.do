cd "x"

******************************************
//Create efficiency variable using merged file
******************************************

use merged_all, clear



******************************************
//Create efficiency variable
******************************************

use merged_all, clear
	

gen max_poss_avg_payoff=(10-(0.5*6+0.5*2))/2  
gen efficiency=(Profit-1.6)/(max_poss_avg_payoff-1.6)   	

sort tgroup markt

save merged_all_eff, replace

collapse efficiency tgroup markt, by(session_id)

save collapsed_all, replace

******************************
//Run Mann-Whitney-tests
******************************

use collapsed_all, clear

ttest efficiency, by(tgroup)
ranksum efficiency, by(tgroup)
