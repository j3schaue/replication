###------------------------------------------------------------###
###------------------------------------------------------------###
###------------------------------------------------------------###
### ANALYSES OF RPP (OSF) REPLICATIONS
###------------------------------------------------------------###
###------------------------------------------------------------###
###------------------------------------------------------------###

# set the working directory to src so we can use relative paths
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(dplyr)
source("../package/replicationTest.R")
source("../package/mdh.R")


###------------------------------------------------------------###
### RPP (OSF)
###------------------------------------------------------------###
data = read.csv("../data/rpp.csv")

experiments = unique(data$experiment)
ks=rep(2, length(experiments))
tau0s = c(0, 1/4, 1/3, 2/3) # plausible ratios for tau0

fedata = lapply(tau0s, FUN=function(tau0)
  setNames(data.frame(
    matrix(unlist(lapply(seq_along(experiments), FUN=function(i) 
      combineResults(t=filter(data, experiment==experiments[i])$z,
                     v=filter(data, experiment==experiments[i])$vz,
                     lambda0=(ks[i]-1)*tau0, 
                     maxratio=20)
    )
    ), ncol=5, byrow = T)
  ), c("k", "Q", paste0("calpha", round(tau0*100, 0)), 
       paste0("p", round(tau0*100, 0)), paste0("mdh", round(tau0*100, 0)))
  )
)


dataout = Reduce(left_join, fedata)
dataout$experiment = experiments
dataout$paper = 'rpp'

dataout = dataout[c("paper", "experiment", "k", "Q", "calpha0",  "p0", "mdh0", "calpha25", "p25",  "mdh25", "calpha33", "p33",  "mdh33", "calpha67", "p67", 
                  "mdh67")] %>%
          left_join(., dplyr::select(data, experiment, replicated, cirep, meta, es, exp_name)) %>%
          distinct()

tmp = data %>% 
  mutate(replicate = as.integer(grepl("_rep", site))) %>%
  group_by(experiment) %>% 
  arrange(replicate) %>%
  summarize(t1 = sum(abs(1 - replicate)*z), t2 = sum(abs(replicate)*z),
            v1=sum(abs(1 - replicate)*vz), v2=sum(abs(replicate)*vz))

dataout = left_join(dataout, tmp)

# write full results
write.csv(dataout, "./results/qtest_fixed_rpp.csv", row.names=F)

# write standard output for combination script
write.csv(dplyr::select(dataout, paper, experiment, k, Q, 
                        calpha0, p0, mdh0, calpha25, p25, mdh25, 
                        calpha33, p33, mdh33, calpha67, p67, mdh67, 
                        t1, t2, v1, v2, replicated),
          "./results/comparison_rpp.csv", row.names=F)

