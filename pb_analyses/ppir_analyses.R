###------------------------------------------------------------###
###------------------------------------------------------------###
###------------------------------------------------------------###
### ANALYSES OF PPIR
###------------------------------------------------------------###
###------------------------------------------------------------###
###------------------------------------------------------------###

# set the working directory to src so we can use relative paths
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(dplyr)
source("../package/replicationTest.R")
source("../package/mdh.R")

# Get data
df = read.csv("../data/ppir.csv")

# ###------------------------------------------------------------###
# ### PPIR Comparison
# ###------------------------------------------------------------###
# #---Comparison Framework--------------------------------
# ## Take the k-1 replicates and create one synthetic estimator
# ## using either fixed or random effects meta-analysis.
# ## Compare the initial finding to this synthetic estimator
# 
# library(metafor)
# 
# # Use a fixed and random effects meta-analysis to combine replications
# methods = c('FE', 'DL')
# 
# # Set the null hypotheses lambda0 = 0, (k-1)/4, (k-1)/3, 2(k-1)/3
# ratios = c(0, 1/4, 1/3, 2/3)
# 
# for(mm in methods){# loop through methods
#   
#   # comp is a table for each lambda0 and each experiment for a given synthetic replicate.
#   comp = lapply(ratios, FUN=function(rr){
#     
#     lambda0 = rr # set the null hypothesis
#     
#     setNames(data.frame(matrix(unlist( # store the results in a data.frame
#       lapply(unique(df$experiment), FUN=function(expt){ # loop through the experiments
#         
#         # separate out the replicates and original study
#         replicates = df %>% filter(experiment==expt & site!='original')
#         orig = df %>% filter(experiment==expt & site=='original')
#         
#         if(nrow(orig) == 1){
#           # combine the replicates
#           tmp = rma.uni(yi=replicates$d, vi=replicates$vd, method=mm)
#           
#           # run a Q-test and get the MDH
#           combineResults(t=c(tmp$beta, orig$d), 
#                          v=c(tmp$se^2, orig$vd),
#                          lambda0=lambda0)
#         } else { list(k=NA, Q=NA, calpha=NA, p=NA, mdh=NA) }
#       })), 
#       ncol=5, byrow=T)), c('k', 'Q', # set the names of the data.frame 'comp'
#                            paste0('calpha', round(100*rr, 0)), 
#                            paste0('p', round(100*rr, 0)), 
#                            paste0('mdh', round(100*rr, 0)))
#     )
#   })
#   
#   # combine all experiments
#   comptab = Reduce(left_join, comp)
#   comptab$experiment = unique(df$experiment)
#   
#   # clean up order of columns and write to file
#   comptab = comptab[c("experiment", "k", "Q", "calpha0",  "p0", "mdh0", "calpha25", "p25",  "mdh25", "calpha33", "p33",  "mdh33", "calpha67", "p67", 
#                       "mdh67")] %>% 
#     left_join(., select(df, experiment, replicated)) %>% distinct()
#   
#   write.csv(comptab, paste0("./results/comparison_ppir_", mm, ".csv"), row.names=F)
# }


###------------------------------------------------------------###
###------------------------------------------------------------###
### Pipeline Project Heterogenetiy Q-Test
###------------------------------------------------------------###
###------------------------------------------------------------###

## Set parameters for analysis
experiments = unique(df$experiment) # unique experiment names
ks = sapply(experiments, # # of trials per experiment
            FUN=function(ee) count(filter(select(df, -n), experiment==ee))$n)
tau0s = c(0, 1/4, 1/3, 2/3) # plausible ratios for tau0
lambda0s = (ks-1)*tau0s  # convert to lambda0
vbars = sapply(experiments, FUN=function(ee) mean(dplyr::filter(df, experiment==ee)$vd)) # avg sampling variances

###----Include original study-----------------------------------
fe = lapply(tau0s, FUN=function(tau0) # loop through null hypotheses lambda0 = 0, (k-1)/4, (k-1)/3, 2(k-1)/3
  setNames(data.frame( # store results as a data frame
    matrix(unlist(lapply(seq_along(experiments), FUN=function(i)
      
      # get the results of the Q-test using all of the studies (rather than aggregating replicates)
      combineResults(t=filter(df, experiment==experiments[i])$d,
                     v=filter(df, experiment==experiments[i])$vd,
                     lambda0=(ks[i]-1)*tau0, 
                     maxratio=100)
    )
    ), ncol=5, byrow = T)
  ), c("k", "Q", paste0("calpha", round(tau0*100, 0)), # name the columns
       paste0("p", round(tau0*100, 0)), paste0("mdh", round(tau0*100, 0)))
  )
)

# Join for all null hypotheses
fetab = Reduce(left_join, fe)
fetab$experiment = experiments
fetab$vdbar = vbars

# reorder df and write to file
fetab = fetab[c("experiment", "k", "Q", "calpha0",  "p0", "mdh0", "calpha25", "p25",  "mdh25", "calpha33", "p33",  "mdh33", "calpha67", "p67", 
                "mdh67", "vbar")]
fetab
write.csv(fetab, "./results/qtest_fixed_ppir_include.csv")


###----Exclude original study-----------------------------------------------
fe = lapply(tau0s, FUN=function(tau0) # loop through null hypotheses lambda0 = 0, (k-1)/4, (k-1)/3, 2(k-1)/3
  setNames(data.frame(  # save as a df
    matrix(unlist(lapply(seq_along(experiments), FUN=function(i) 
      
      # get Q-test and MDH excluding the original study
      combineResults(t=filter(df, experiment==experiments[i] & site!='original')$d,
                     v=filter(df, experiment==experiments[i] & site!='original')$vd,
                     lambda0=(ks[i]-1)*tau0, 
                     maxratio=100)
    )
    ), ncol=5, byrow = T)
  ), c("k", "Q", paste0("calpha", round(tau0*100, 0)), 
       paste0("p", round(tau0*100, 0)), paste0("mdh", round(tau0*100, 0)))
  )
)

# Join for all null hypotheses
fetab = Reduce(left_join, fe)
fetab$experiment = experiments
fetab$vdbar = vbars

# reorder df and write to file
fetab = fetab[c("experiment", "k", "Q", "calpha0",  "p0", "mdh0", "calpha25", "p25",  "mdh25", "calpha33", "p33",  "mdh33", "calpha67", "p67", 
                "mdh67", "vbar")]
fetab
write.csv(fetab, "./results/qtest_fixed_ppir_exclude.csv")


###------------------------------------------------------------###
###------------------------------------------------------------###
### Pipeline Project Variance Components
###------------------------------------------------------------###
###------------------------------------------------------------###

# Estimate variance components with Paule-Mandel and DerSimonian-Laird estimators
methods = c("PM", "DL")

###----Include Original Study--------------------------------------------------
for(mm in methods){# for each method, compute tau^2 for each set of replicates
  tab = as.data.frame(do.call(rbind, lapply(experiments, FUN=function(ee){
    dd = filter(df, experiment==ee)
    ff = confint(rma.uni(yi=dd$d, vi=dd$vd, method=mm))$random[1,] # get the estimate and the CI
  })))
  tab$experiment = experiments
  tab$vdbar = vbars
  write.csv(select(tab, experiment, tau2=estimate, ci.lb, ci.ub), paste0('./results/ppir_vc_include_', mm,'.csv'), row.names=F)
}

###----Exclude Original Study--------------------------------------------------
for(mm in methods){# for each method, compute tau^2 for each set of replicates
  tab = as.data.frame(do.call(rbind, lapply(experiments, FUN=function(ee){
    dd = filter(df, experiment==ee & site!='original') # exclude the original study
    ff = confint(rma.uni(yi=dd$d, vi=dd$vd, method=mm))$random[1,] # get the estimate and the CI
  })))
  tab$experiment = experiments
  tab$vdbar = vbars
  write.csv(select(tab, experiment, tau2=estimate, ci.lb, ci.ub), paste0('./results/ppir_vc_exclude_', mm,'.csv'), row.names=F)
}