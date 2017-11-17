###------------------------------------------------------------###
###------------------------------------------------------------###
###------------------------------------------------------------###
### ANALYSES OF Alogna (RRR)
###------------------------------------------------------------###
###------------------------------------------------------------###
###------------------------------------------------------------###

# set the working directory to src so we can use relative paths
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(dplyr)
source("../package/replicationTest.R")
source("../package/mdh.R")

# Get data
df = read.csv("../data/rrr_alogna.csv")
experiments = unique(data$experiment)

###------------------------------------------------------------###
### Alogna Comparison
###------------------------------------------------------------###

#---Fixed effects meta-analyses
library(metafor)
ratios = c(0, 1/4, 1/3, 2/3)
methods = c('FE', 'DL')
for(mm in methods){
  comp = lapply(ratios, FUN=function(rr){
    lambda0 = rr
    setNames(data.frame(matrix(unlist(
      lapply(unique(df$experiment), FUN=function(expt){
        replicates = df %>% filter(experiment==expt & site!='original' & !is.infinite(rd))
        orig = df %>% filter(experiment==expt & site=='original' & !is.na(rd))
        if(nrow(orig) == 1){
          tmp = rma.uni(yi=replicates$rd, vi=replicates$vrd, method='FE')
          combineResults(t=c(tmp$beta, orig$rd), 
                         v=c(tmp$se^2, orig$vrd),
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
  comptab = comptab[c("experiment", "k", "Q", "calpha0",  "p0", "mdh0", "calpha25", "p25",
                      "mdh25", "calpha33", "p33",  "mdh33", "calpha67", "p67", "mdh67")] %>% 
    left_join(., select(df, experiment, replicated)) %>% distinct()
  
  write.csv(comptab, paste0("./results/comparison_alogna_", mm, ".csv"), row.names=F)
}
