###------------------------------------------------------------###
###------------------------------------------------------------###
###------------------------------------------------------------###
### ANALYSES OF RPP (OSF) REPLICATIONS
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
### RPP (OSF)
###------------------------------------------------------------###
data = read.csv("../data/rpp.csv")

experiments = unique(data$experiment)
# for(ee in experiments){
#   dd = filter(data, experiment==ee)
#   n = nrow(dd)
#   show(replicationTest(t=dd$z, v=dd$vz))
#   show(mdh_constvar(k=n, maxratio = 100))
# }

ks=rep(2, length(experiments))

fedata = lapply(tau0s, FUN=function(tau0)
  setNames(data.frame(
    matrix(unlist(lapply(seq_along(experiments), FUN=function(i) 
      combineResults(t=filter(data, experiment==experiments[i])$z,
                     v=filter(data, experiment==experiments[i])$vz,
                     lambda0=(ks[i]-1)*tau0, 
                     maxratio=100)
    )
    ), ncol=5, byrow = T)
  ), c("k", "Q", paste0("calpha", round(tau0*100, 0)), 
       paste0("p", round(tau0*100, 0)), paste0("mdh", round(tau0*100, 0)))
  )
)


dataout = Reduce(left_join, fedata)
dataout$experiment = experiments
dataout = dataout[c("experiment", "k", "Q", "calpha0",  "p0", "mdh0", "calpha25", "p25",  "mdh25", "calpha33", "p33",  "mdh33", "calpha67", "p67", 
                  "mdh67")]

write.csv(dataout, "./results/qtest_fixed_rpp.csv", row.names=F)
round(camout$p0, 3)