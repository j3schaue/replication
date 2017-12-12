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
source("../package/mdh.R"); source("./misc.R")
paper = 'rpe'
tes = 'z'; vr = 'vz'


###------------------------------------------------------------###
### RPE (Camerer) Comparison Analyses on z-Scale
###------------------------------------------------------------###
cam = read.csv("../data/rpe.csv")

# re-code sites for analyses
cam$site[grepl("orig", cam$site)] = 'original'

experiments = unique(cam$experiment)
ks=rep(2, length(experiments))
tau0s = c(0, 1/4, 1/3, 2/3)

# Use a fixed and random effects meta-analysis to combine replications
methods = c('FE', 'DL')

# Set the null hypotheses lambda0 = 0, (k-1)/4, (k-1)/3, 2(k-1)/3
ratios = c(0, 1/4, 1/3, 2/3)

fecam = lapply(ratios, FUN=function(tau0)
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
  left_join(., dplyr::select(cam, experiment, replicated, es)) %>%
  distinct()

tmp = cam %>% group_by(experiment) %>% arrange(replicate) %>%
    summarize(t1 = sum(abs(1 - replicate)*z), t2 = sum(abs(replicate)*z),
              v1=sum(abs(1 - replicate)*vz), v2=sum(abs(replicate)*vz))
camout = left_join(camout, tmp)

camout$es = 'z'
camout$paper = 'rpe'

camout

# write standard output for combination script
write.csv(dplyr::select(camout, paper, experiment, k, Q, 
                        calpha0, p0, mdh0, calpha25, p25, mdh25, 
                        calpha33, p33, mdh33, calpha67, p67, mdh67, 
                        t1, t2, v1, v2, replicated),
          "./results/comparison_rpe_FE.csv", row.names=F)

###------------------------------------------------------------###
### RPE (Camerer) Comparison Analyses on Cohen's d Scale
###------------------------------------------------------------###

dcam = lapply(tau0s, FUN=function(tau0)
  setNames(data.frame(
    matrix(unlist(lapply(seq_along(experiments), FUN=function(i) 
      combineResults(t=filter(cam, experiment==experiments[i])$d,
                     v=filter(cam, experiment==experiments[i])$vd,
                     lambda0=(ks[i]-1)*tau0,
                     maxratio=20)
    )
    ), ncol=5, byrow = T)
  ), c("k", "Q", paste0("calpha", round(tau0*100, 0)), 
       paste0("p", round(tau0*100, 0)), paste0("mdh", round(tau0*100, 0)))
  )
)

dout = Reduce(left_join, dcam)
dout$experiment = experiments
dout = dout[c("experiment", "k", "Q", "calpha0",  "p0", "mdh0", "calpha25", "p25",  "mdh25", "calpha33", "p33",  "mdh33", "calpha67", "p67", 
              "mdh67")] %>%
  left_join(., dplyr::select(cam, experiment, replicated, es)) %>%
  distinct()

tmp = cam %>% group_by(experiment) %>% arrange(replicate) %>%
  summarize(t1 = sum(abs(1 - replicate)*d), t2 = sum(abs(replicate)*d),
            v1=sum(abs(1 - replicate)*vd), v2=sum(abs(replicate)*vd))
dout = left_join(dout, tmp)

dout$es = 'd'
dout$paper = 'rpe'

# write standard output for combination script
write.csv(dplyr::select(dout, paper, experiment, k, Q, 
                        calpha0, p0, mdh0, calpha25, p25, mdh25, 
                        calpha33, p33, mdh33, calpha67, p67, mdh67, 
                        t1, t2, v1, v2, replicated),
          "./results/comparison_rpe_FE_d.csv", row.names=F)
