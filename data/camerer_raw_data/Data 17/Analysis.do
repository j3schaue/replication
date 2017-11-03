clear
set more off

use "Data.dta"

* Drop Practice Round
drop if period == 0

sort    market period subject
rename  os2effort output

* Generate Unique Identifier
egen uid = group(session group)
collapse output market groupsize insider group, by(uid period)

* Group output: Insiders
gllamm output market period if insider==1, i(uid) link(oprobit) robust
test _b[market]=0
