###-----------------------------------------------------------###
###-----------------------------------------------------------###
### Table 1S from Camerer
###  Reports effect sizes (r) and sample sizes (n) per lab
###  Cobbled from PDF in paper.
###-----------------------------------------------------------###
###-----------------------------------------------------------###

setwd("~/Documents/replication/replication/data/camerer_raw_data/")
list.files()
source("S1pvalsetc.R")
orig = read.csv("S1origes.csv")
reps = read.csv("S1repes.csv")
exps = read.csv("S1experiment.csv")
repdet = read.csv("S1replicationdet.csv")
s1 = do.call(cbind, list(exps, orig, reps, repdet))

write.csv(s1, "../camerer_S1.csv", row.names=F)
