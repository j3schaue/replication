###------------------------------------------------------------###
###------------------------------------------------------------###
###------------------------------------------------------------###
### ANALYSES OF RPE REPLICATIONS
###------------------------------------------------------------###
###------------------------------------------------------------###
###------------------------------------------------------------###

# set the working directory to src so we can use relative paths
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(dplyr)
source("../package/replicationTest.R")
source("../package/mdh.R")


###------------------------------------------------------------###
### RPE (Camerer)
###------------------------------------------------------------###
cam = read.csv("../data/rpe.csv")

experiments = unique(cam$experiment)
ks=rep(2, length(experiments))
tau0s = c(0, 1/4, 1/3, 2/3)

fecam = lapply(tau0s, FUN=function(tau0)
  setNames(data.frame(
    matrix(unlist(lapply(seq_along(experiments), FUN=function(i) 
      combineResults(t=filter(cam, experiment==experiments[i])$z,
                     v=filter(cam, experiment==experiments[i])$vz,
                     lambda0=(ks[i]-1)*tau0,
                     maxratio=20)
    )
    ), ncol=5, byrow = T)
  ), c("k", "Q", paste0("calpha", round(tau0*100, 0)), 
       paste0("p", round(tau0*100, 0)), paste0("mdh", round(tau0*100, 0)))
  )
)

camout = Reduce(left_join, fecam)
camout$experiment = experiments
camout = camout[c("experiment", "k", "Q", "calpha0",  "p0", "mdh0", "calpha25", "p25",  "mdh25", "calpha33", "p33",  "mdh33", "calpha67", "p67", 
                  "mdh67")] %>%
  left_join(., dplyr::select(cam, experiment, replicated, es, vbar, v0)) %>%
  distinct()
camout$es = 'z'
camout$paper = 'rpe'
tmp = cam %>% group_by(experiment) %>% arrange(replicate) %>%
        summarize(vbar=mean(vz), v0=sum(abs(1 - replicate)*vz))
camout = left_join(camout, tmp)

# write full results
write.csv(camout, "./results/qtest_fixed_rpe.csv", row.names=F)

# write standard output for combination script
write.csv(dplyr::select(camout, paper, experiment, k, Q, 
                        calpha0, p0, mdh0, calpha25, p25, mdh25, 
                        calpha33, p33, mdh33, calpha67, p67, mdh67, 
                        vbar, replicated),
          "./results/qtest_fixed_rpe_include.csv", row.names=F)
