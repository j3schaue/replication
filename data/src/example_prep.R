###--------------------------------------------###
###--------------------------------------------###
### DATA PREPARATION FOR EXAMPLES
###  This code preps the data from some ManyLabs
###   replications in order to use them as examples
###   for the replication project.
###
### We use only the Gambler's Fallacy example.
### Collin will add other data.
###--------------------------------------------###
###--------------------------------------------###


###--------------------------------------------###
### GAMBLER'S FALLACY
###--------------------------------------------###

# read in the data and remove the 'notes' columns
library(xlsx)
# gfa <- read.xlsx("ml_gamblers-fallacy.xlsx", sheetIndex=1, header=T)
gfa = read.xlsx("ML-_Summary_Statistics.xlsx", "Gambler's Fallacy")
dim(gfa)
gfa <- gfa[, 1:(ncol(gfa) - 2)] # drop the summary columns
names(gfa)

# give the columns easy to read names
names(gfa) = c("Site", "NThree6", "NTwo6", "NExcluded", "MeanThree6", "MeanTwo6", "SDThree6", "SDTwo6", "tev", "tuev", "dfev", "dfuev", "ESuev", "ESmd")
head(gfa)

# store the overall results in a separate df
gfa_overall <- gfa[1:2, ]
gfa <- gfa[3:nrow(gfa), ]
str(gfa)


###--------------------------------------------###
### Check to see how the ES and other columns 
### are calculated
###--------------------------------------------###
# Check the equal variance for the t-statistics
pooled_var <- ((gfa$NThree6 - 1) * gfa$SDThree6^2 + (gfa$NTwo6 - 1) * gfa$SDTwo6^2)/(gfa$NThree6 + gfa$NTwo6 - 2) # pooled variance
mytev <- (gfa$MeanThree6 - gfa$MeanTwo6)/(sqrt(pooled_var) * sqrt(1/gfa$NThree6 + 1/gfa$NTwo6)) # my calculated t-stats
sort((gfa$tev - mytev)/gfa$tev) * 100 # percent difference between their and mine are all about 0

# Check the unequal variance t-statistics
mytuev <- (gfa$MeanThree6 - gfa$MeanTwo6)/sqrt(gfa$SDThree6^2/gfa$NThree6 + gfa$SDTwo6^2/gfa$NTwo6) # my calculated t-stats
sort(gfa$tuev - mytuev) * 100 # percent difference between their and mine are all about 0

# Check their effect sizes
bad_pooled_var <- (gfa$SDThree6^2 + gfa$SDTwo6^2)/2 # It looks like they just took the mean variance for the pooled variance estimate
myesmd_bad <- (gfa$MeanThree6 - gfa$MeanTwo6)/sqrt(bad_pooled_var) # My (incorrect) SMDs
(gfa$ESmd - myesmd_bad)/gfa$ESmd * 100 # percent difference indicates that they use just a mean variance rather than a pooled variance


###------------------------------------------###
### Compute the effect sizes we want
###------------------------------------------###
gfa$J <- 1 - 3/(4*gfa$dfev - 1) # small sample correction for g
gfa$d <- (gfa$MeanThree6 - gfa$MeanTwo6)/sqrt(pooled_var) # SMD, Cohen's d, pooled estimate of common variance calculated above
gfa$g <- gfa$d * gfa$J # Glass's (Hedges') g

gfa$vd <- (gfa$NThree6 + gfa$NTwo6) / (gfa$NThree6 * gfa$NTwo6) + gfa$d^2/(2 * (gfa$NThree6 + gfa$NTwo6)) # Variance of d
gfa$vg <- gfa$vd * gfa$J^2 # Variance of g

# Add original study
toadd = data.frame(Site='Original', d=0.69, vd=0.0718, dfev=57, NThree6=30, NTwo6=29)
toadd$J <- 1 - 3/(4*toadd$dfev - 1)
toadd$g <- toadd$d * toadd$J
toadd$vg <- toadd$vd * toadd$J^2 # Variance of g
gfa = dplyr::bind_rows(toadd, gfa)
write.csv(gfa, "gf_example_fixed.csv", row.names=F)

source("replicationTest.R")
qStat(t=gfa$d, v=gfa$vd)
replicationTest(t=gfa$d, v=gfa$vd, burden_on_replication = F, fixed=T, exact = F, lambda0=35/3)$p
