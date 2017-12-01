clear
set more off

use "Data.dta"

collapse (mean) profitprincipal profitagent, by(treatment session type)
gen     surplus = profitprincipal + profitagent

mean    profitprincipal if type==1 & treatment==0
mean    profitprincipal if type==1 & treatment==1 
mean    profitagent     if type==2 & treatment==0 
mean    profitagent     if type==2 & treatment==1
mean    surplus         if type==1 & treatment==0
mean    surplus         if type==1 & treatment==1

ranksum profitprincipal if type==1, by(treatment) 
ranksum profitagent     if type==2, by(treatment)
ranksum surplus         if type==1, by(treatment) 

graph bar profitprincipal profitagent surplus, by(treatment)
