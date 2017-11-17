###------------------------------------------------------------###
###------------------------------------------------------------###
###------------------------------------------------------------###
### ANALYSES OF MANYLABS REPLICATIONS
###------------------------------------------------------------###
###------------------------------------------------------------###
###------------------------------------------------------------###

# set the working directory to src so we can use relative paths
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(dplyr)
source("../package/replicationTest.R")
source("../package/mdh.R")

###------------------------------------------------------------###
### Many Labs Comparison
###------------------------------------------------------------###
df = read.csv("../data/manylabs_comparison.csv")

#---Fixed effects meta-analyses
library(metafor)
ratios = c(0, 1/4, 1/3, 2/3)
methods = c('FE', 'DL')
for(mm in methods){
  comp = lapply(ratios, FUN=function(rr){
    lambda0 = rr
    setNames(data.frame(matrix(unlist(
      lapply(unique(df$experiment), FUN=function(expt){
      replicates = df %>% filter(experiment==expt & site!='original' & !is.infinite(t))
      orig = df %>% filter(experiment==expt & site=='original' & !is.na(t))
      if(nrow(orig) == 1){
        tmp = rma.uni(yi=replicates$t, vi=replicates$v, method='FE')
        combineResults(t=c(tmp$beta, orig$t), 
                     v=c(tmp$se^2, orig$v),
                     lambda0=lambda0)
      } else { list(k=NA, Q=NA, calpha=NA, p=NA, mdh=NA) }
    })), 
    ncol=5, byrow=T)), c('k', 'Q', 
                       paste0('calpha', round(100*rr, 0)), 
                       paste0('p', round(100*rr, 0)), 
                       paste0('mdh', round(100*rr, 0)))
    )
  })
  
  comptab = Reduce(left_join, comp)
  comptab$experiment = unique(df$experiment)
  comptab = comptab[c("experiment", "k", "Q", "calpha0",  "p0", "mdh0", "calpha25", "p25",  "mdh25", "calpha33", "p33",  "mdh33", "calpha67", "p67", 
                  "mdh67")] %>% 
            left_join(., select(df, experiment, replicated)) %>% distinct()
  
  write.csv(comptab, paste0("./results/comparison_manylabs_", mm, ".csv"), row.names=F)
}

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
