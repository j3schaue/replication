###------------------------------------------------------------###
###------------------------------------------------------------###
###------------------------------------------------------------###
### ANALYSES OF MULTI-LAB REPLICATION EFFORTS
###------------------------------------------------------------###
###------------------------------------------------------------###
###------------------------------------------------------------###

library(dplyr)
setwd("~/Documents/replication/replication/pb_analyses/")
source("../package/replicationTest.R")
source("../package/mdh.R")
combineResults = function(t=NULL, v=NULL, h0replication=TRUE, fixed=TRUE, alpha=.05, lambda0=0, tau0=0, power=0.8, step=.001, maxratio=100){
  qtest = replicationTest(t=t, v=v, h0replication=h0replication, fixed=fixed, alpha=alpha, lambda0=lambda0, tau0=tau0)
  qtest[["mdh"]] = mdh_constvar(k=length(t), alpha=alpha, power=power, h0replication=h0replication, lambda0=lambda0, step=step, maxratio=maxratio)
  return(qtest[c("k", "Q", "calpha", "p", "mdh")])
}

###------------------------------------------------------------###
### Many Labs
###------------------------------------------------------------###
## Read in data
data = read.csv("../data/Data_Clean.csv") %>% filter(is.finite(g)) # drop infinite estimates: Collin--are these just missing values?
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
fetab
# fetabl$mdh_ratio = sapply(ks, FUN=function(k) mdh_constvar(k))
# f_ml$mdh = f_ml$mdh_ratio * vbars