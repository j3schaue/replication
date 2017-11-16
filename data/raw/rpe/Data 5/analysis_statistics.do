use "20150824_data_effort.dta", clear

generate ingroup = 0
generate outgroup = 0
replace ingroup = 1 if treatment == 1
replace outgroup = 1 if treatment == -1
generate inenh = ingroup*enhanced
generate outenh = outgroup*enhanced

generate age18 = 0
replace age18 = 1 if age == 18
generate age19 = 0
replace age19 = 1 if age == 19
generate age20 = 0
replace age20 = 1 if age == 20
generate age21 = 0
replace age21 = 1 if age == 21
generate age22 = 0
replace age22 = 1 if age == 22

generate raceAsian = 0
replace raceAsian = 1 if race == 0
generate raceBlack = 0
replace raceBlack = 1 if race == 1
generate raceHispanic = 0
replace raceHispanic = 1 if race == 3
generate raceMulti = 0
replace raceMulti = 1 if race == 5
generate raceOther = 0
replace raceOther = 1 if race == 6

generate raisedAsia = 0
replace raisedAsia = 1 if country == 2

generate married = 0
replace married = 1 if marriage == 1

generate employedFull = 0
replace employedFull = 1 if employment == 0
generate employedPart = 0
replace employedPart = 1 if employment == 1

generate onesib = 0
replace onesib = 1 if siblings == 1
generate twosib = 0
replace twosib = 1 if siblings == 2
generate moresib = 0
replace moresib = 1 if siblings >= 3

*Effort regression: group salience (Table 2)
*No Demographics
xtreg effort inenh outenh, re i(subject) cluster(session) 
test inenh == outenh

*Demographics
xtreg effort inenh outenh age18 age19 age20 age21 age22 gender raceAsian raceBlack raceHispanic raceMulti raceOther raisedAsia married employedFull employedPart onesib twosib moresib expensessharedspouse-expensesother voted-volunteer, re i(subject) cluster(session)
test inenh == outenh
