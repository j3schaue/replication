clear
set more off

use "Data.dta"

collapse (mean) m6 (max) period (sum) earn_per, by(session)

* Calculate Efficiency Ratio
local v_eff = (7 * log(7) - 6)
gen   v_eff = period * 3 * `v_eff'
gen   eff_ratio = earn_per / v_eff

sum   eff_ratio if m6 == 0
sum   eff_ratio if m6 == 1

ranksum eff_ratio, by(m6)
