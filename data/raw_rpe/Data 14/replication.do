/* Replication of Ifcher and Zarghamee (2011) AER */

xi: reg pv vh1 vh3 i.fv*i.delay i.college i.female i.race i.religion i.practicing i.family_income i.happiness if pv<fv & treatment == 0, cluster(id)

