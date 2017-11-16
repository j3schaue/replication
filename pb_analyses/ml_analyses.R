###------------------------------------------------------------###
###------------------------------------------------------------###
###------------------------------------------------------------###
### ANALYSES OF MANYLABS REPLICATIONS
###------------------------------------------------------------###
###------------------------------------------------------------###
###------------------------------------------------------------###

library(dplyr)
setwd("~/Documents/replication/replication/pb_analyses/") #jakes path
# setwd("./pb_analyses") # relative path
source("../package/replicationTest.R")
source("../package/mdh.R")
combineResults = function(t=NULL, v=NULL, h0replication=TRUE, fixed=TRUE, alpha=.05, lambda0=0, tau0=0, power=0.8, step=.001, maxratio=100){
  qtest = replicationTest(t=t, v=v, h0replication=h0replication, fixed=fixed, alpha=alpha, lambda0=lambda0, tau0=tau0)
  qtest[["mdh"]] = mdh_constvar(k=length(t), alpha=alpha, power=power, h0replication=h0replication, lambda0=lambda0, step=step, maxratio=maxratio)
  return(qtest[c("k", "Q", "calpha", "p", "mdh")])
}

###------------------------------------------------------------###
### Many Labs Comparison
###------------------------------------------------------------###




###------------------------------------------------------------###
### Many Labs pregistered
###------------------------------------------------------------###
## Read in data
data = read.csv("../data/manylabs.csv") %>% filter(is.finite(g)) # drop infinite estimates: Collin--are these just missing values?
names(data)[names(data) == "es.measurement"] = "es" # simplify names

## Set parameters for analysis
experiments = unique(data$Experiment) # unique experiment names
ks = sapply(experiments, # # of trials per experiment
            FUN=function(ee) count(filter(data, Experiment==ee))$n)
tau0s = c(0, 1/4, 1/3, 2/3) # plausible ratios for tau0
lambda0s = (ks-1)*tau0s  # convert to lambda0
vbars = sapply(experiments, FUN=function(ee) mean(filter(data, Experiment==ee)$vg)) # avg sampling variances


fe = lapply(tau0s, FUN=function(tau0)
  setNames(data.frame(
    matrix(unlist(lapply(seq_along(experiments), FUN=function(i) 
      combineResults(t=filter(data, Experiment==experiments[i])$g,
                  v=filter(data, Experiment==experiments[i])$vg,
                  lambda0=(ks[i]-1)*tau0, 
                  maxratio=100)
      )
    ), ncol=5, byrow = T)
  ), c("k", "Q", paste0("calpha", round(tau0*100, 0)), 
       paste0("p", round(tau0*100, 0)), paste0("mdh", round(tau0*100, 0)))
  )
)

fetab = Reduce(left_join, fe)
fetab$experiment = experiments
fetab = fetab[c("experiment", "k", "Q", "calpha0",  "p0", "mdh0", "calpha25", "p25",  "mdh25", "calpha33", "p33",  "mdh33", "calpha67", "p67", 
                "mdh67")]
write.csv(fetab, "./results/qtest_fixed_manylabs.csv")


# ###--- Random Effects Analyses
# redata = lapply(tau0s, FUN=function(tau0)
#   setNames(data.frame(
#     matrix(unlist(lapply(seq_along(experiments), FUN=function(i) 
#       combineResults(t=filter(data, experiment==experiments[i])$d,
#                      v=filter(data, experiment==experiments[i])$vd,
#                      lambda0=(ks[i]-1)*tau0, 
#                      fixed=FALSE,
#                      maxratio=100)
#     )
#     ), ncol=5, byrow = T)
#   ), c("k", "Q", paste0("calpha", round(tau0*100, 0)), 
#        paste0("p", round(tau0*100, 0)), paste0("mdh", round(tau0*100, 0)))
#   )
# )
# 
# 
# rdataout = Reduce(left_join, redata)
# rdataout$experiment = experiments
# rdataout = rdataout[c("experiment", "k", "Q", "calpha0",  "p0", "mdh0", "calpha25", "p25",  "mdh25", "calpha33", "p33",  "mdh33", "calpha67", "p67", 
#                       "mdh67")]
# 
# write.csv(rdataout, "./results/qtest_random_ppir.csv", row.names=F)
# round(rdataout$p0, 3)
