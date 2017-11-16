## Computation of original effect sizes
# The manylabs.xlsx workbook contains info on original effect sizes
# reported as Cohen's d. However, this may not be accurate, nor are
# there sampling variances provided.
# This script documents efforts to compute sampling variances.
# In doing so, it became evident that some original effect sizes 
# ought to be recomputed.

#---------------------------------------#
# Effect size for Tversky from p. 1
# Gain-Loss
#---------------------------------------#
# data
n1 = 152; n2 = 155
n11 = n1*.72; n12 = n1*.28 # Problem 1
n21 = n2*.22; n22 = n2*.78 # Problem 2
data = round(c(n11, n12, n21, n22), 0)

# log odds ratio
or = n11*n22/(n21*n12)
log_or = log(or)
vlog_or = 1/n11 + 1/n22 + 1/n21 + 1/n12
d = sqrt(3)/pi * log_or
vd = 3/pi * vlog_or

#---------------------------------------#
# Effect size for Schwartz from Table 1
# Scales
#---------------------------------------#
# data
n1 = 68; n2 = 64
n11 = n1*(1-.162); n12 = n1*(.162) # Problem 1
n21 = n2*(.625); n22 = n2*(1-.625) # Problem 2
data = round(c(n1*c(1-.162, .162), n2*c(.625, 1-.625)), 0)

# log odds ratio
or = n11*n22/(n21*n12)
log_or = log(or)
vlog_or = 1/n11 + 1/n22 + 1/n21 + 1/n12
d = sqrt(3)/pi * log_or
vd = 3/pi * vlog_or

#---------------------------------------#
# Effect size for Carter from p. 1061
#---------------------------------------#
n = 66; nt = 33; nc = 33
mt = 2.65; mc = 3.10 
tstat = -2.04
tstat / sqrt(nt*nc/(nt + nc))


#---------------------------------------#
# Effect size for Oppenheimer 2009 
# from p. 1061
#---------------------------------------#
n = 213
mt = 7.46; mc = 6.93 
fstat = 2.74
tstat = sqrt(2.74)
tstat / sqrt(nt*nc/(nt + nc))
tstat/sqrt(n/4)

#---------------------------------------#
# Check for Nosek
#---------------------------------------#
n = 213
rr = .42
z = 0.5 * log((1 + rr)/(1 - rr))
vz = 1/(n-3)

#---------------------------------------#
# Effect size for Rugg 1941
#---------------------------------------#
p11 = .62; p12 = .21; p21 = .39; p22 = .46
log(p11*p22/(p12*p21)) * sqrt(3)/pi
