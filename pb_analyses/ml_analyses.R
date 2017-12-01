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
###------------------------------------------------------------###
### Many Labs Comparison
###------------------------------------------------------------###
###------------------------------------------------------------###
# load the data
df = read.csv("../data/manylabs_comp.csv") %>% filter(is.finite(t)) # weed out infinite estimates

#---Comparison Framework--------------------------------
## Take the k-1 replicates and create one synthetic estimator
## using either fixed or random effects meta-analysis.
## Compare the initial finding to this synthetic estimator

library(metafor)

# Use a fixed and random effects meta-analysis to combine replications
methods = c('FE', 'DL')

# Set the null hypotheses lambda0 = 0, (k-1)/4, (k-1)/3, 2(k-1)/3
ratios = c(0, 1/4, 1/3, 2/3)

for(mm in methods){# loop through methods
  
  # comp is a table for each lambda0 and each experiment for a given synthetic replicate.
  comp = lapply(ratios, FUN=function(rr){
    
    lambda0 = rr # set the null hypothesis
    
    setNames(data.frame(matrix(unlist( # store the results in a data.frame
      lapply(unique(df$experiment), FUN=function(expt){ # loop through the experiments
      
      # separate out the replicates and original study
      replicates = df %>% filter(experiment==expt & site!='original')
      orig = df %>% filter(experiment==expt & site=='original')
      
      if(nrow(orig) == 1){
        # combine the replicates
        tmp = rma.uni(yi=replicates$t, vi=replicates$v, method=mm)
        
        # run a Q-test and get the MDH
        combineResults(t=c(tmp$beta, orig$t), 
                     v=c(tmp$se^2, orig$v),
                     lambda0=lambda0, maxratio=20)
      } else { list(k=NA, Q=NA, calpha=NA, p=NA, mdh=NA) }
    })), 
    ncol=5, byrow=T)), c('k', 'Q', # set the names of the data.frame 'comp'
                       paste0('calpha', round(100*rr, 0)), 
                       paste0('p', round(100*rr, 0)), 
                       paste0('mdh', round(100*rr, 0)))
    )
  })
  
  # combine all experiments
  comptab = Reduce(left_join, comp)
  comptab$experiment = unique(df$experiment)
  
  # clean up order of columns and write to file
  comptab = comptab[c("experiment", "k", "Q", "calpha0",  "p0", "mdh0", "calpha25", "p25",  "mdh25", "calpha33", "p33",  "mdh33", "calpha67", "p67", 
                  "mdh67")] %>% 
            left_join(., dplyr::select(df, experiment, replicated)) %>% distinct()
  comptab$paper = 'manylabs'
  
  write.csv(comptab, paste0("./results/comparison_manylabs_", mm, ".csv"), row.names=F)
}


###------------------------------------------------------------###
###------------------------------------------------------------###
### Many Labs Heterogenetiy Q-Test
###------------------------------------------------------------###
###------------------------------------------------------------###

## Set parameters for analysis
experiments = unique(df$experiment) # unique experiment names
ks = sapply(experiments, # no. of trials per experiment
            FUN=function(ee) count(filter(df, experiment==ee))$n)
tau0s = c(0, 1/4, 1/3, 2/3) # plausible ratios for tau0
lambda0s = (ks-1)*tau0s  # convert to lambda0
vbars = sapply(experiments, FUN=function(ee) mean(dplyr::filter(df, experiment==ee)$v)) # avg sampling variances

###----Include original study-----------------------------------
fe = lapply(tau0s, FUN=function(tau0) # loop through null hypotheses lambda0 = 0, (k-1)/4, (k-1)/3, 2(k-1)/3
  setNames(data.frame( # store results as a data frame
    matrix(unlist(lapply(seq_along(experiments), FUN=function(i)
      
      # get the results of the Q-test using all of the studies (rather than aggregating replicates)
      combineResults(t=filter(df, experiment==experiments[i])$t,
                  v=filter(df, experiment==experiments[i])$v,
                  lambda0=(ks[i]-1)*tau0, 
                  maxratio=20)
      )
    ), ncol=5, byrow = T)
  ), c("k", "Q", paste0("calpha", round(tau0*100, 0)), # name the columns
       paste0("p", round(tau0*100, 0)), paste0("mdh", round(tau0*100, 0)))
  )
)

# Join for all null hypotheses
fetab = Reduce(left_join, fe)
fetab$experiment = experiments
fetab$vbar = vbars
fetab$paper = 'manylabs'

# reorder df and write to file
fetab = fetab[c("paper", "experiment", "k", "Q", "calpha0",  "p0", "mdh0", "calpha25", "p25",  "mdh25", "calpha33", "p33",  "mdh33", "calpha67", "p67", 
                "mdh67", "vbar")] %>%
          left_join(., distinct(dplyr::select(df, experiment, replicated)))
fetab
write.csv(fetab, "./results/qtest_fixed_manylabs_include.csv", row.names=F)


###----Exclude original study-----------------------------------------------
fe = lapply(tau0s, FUN=function(tau0) # loop through null hypotheses lambda0 = 0, (k-1)/4, (k-1)/3, 2(k-1)/3
  setNames(data.frame(  # save as a df
    matrix(unlist(lapply(seq_along(experiments), FUN=function(i) 
      
      # get Q-test and MDH excluding the original study
      combineResults(t=filter(df, experiment==experiments[i] & site!='original')$t,
                     v=filter(df, experiment==experiments[i] & site!='original')$v,
                     lambda0=(ks[i]-1)*tau0, 
                     maxratio=20)
    )
    ), ncol=5, byrow = T)
  ), c("k", "Q", paste0("calpha", round(tau0*100, 0)), 
       paste0("p", round(tau0*100, 0)), paste0("mdh", round(tau0*100, 0)))
  )
)

# Join for all null hypotheses
fetab = Reduce(left_join, fe)
fetab$experiment = experiments
fetab$vbar = vbars
fetab$paper = "manylabs"

# reorder df and write to file
fetab = fetab[c("paper", "experiment", "k", "Q", "calpha0",  "p0", "mdh0", "calpha25", "p25",  "mdh25", "calpha33", "p33",  "mdh33", "calpha67", "p67", 
                "mdh67", "vbar")] %>%
  left_join(., distinct(dplyr::select(df, experiment, replicated))) # add replication designation

fetab
write.csv(fetab, "./results/qtest_fixed_manylabs_exclude.csv", row.names=F)


###------------------------------------------------------------###
###------------------------------------------------------------###
### Many Labs variance components
###------------------------------------------------------------###
###------------------------------------------------------------###

# Estimate variance components with Paule-Mandel and DerSimonian-Laird estimators
methods = c("PM", "DL")

###----Include Original Study--------------------------------------------------
for(mm in methods){# for each method, compute tau^2 for each set of replicates
  tab = as.data.frame(do.call(rbind, lapply(experiments, FUN=function(ee){
    dd = filter(df, experiment==ee)
    ff = confint(rma.uni(yi=dd$t, vi=dd$v, method=mm))$random[1,] # get the estimate and the CI
  })))
  tab$experiment = experiments
  tab$vbar = vbars
  tab$paper = 'manylabs'
  write.csv(dplyr::select(tab, experiment, tau2=estimate, ci.lb, ci.ub, paper), paste0('./results/manylabs_vc_include_', mm,'.csv'), row.names=F)
}

###----Exclude Original Study--------------------------------------------------
for(mm in methods){# for each method, compute tau^2 for each set of replicates
  tab = as.data.frame(do.call(rbind, lapply(experiments, FUN=function(ee){
    dd = filter(df, experiment==ee & site!='original') # exclude the original study
    ff = confint(rma.uni(yi=dd$t, vi=dd$v, method=mm))$random[1,] # get the estimate and the CI
  })))
  tab$experiment = experiments
  tab$vbar = vbars
  tab$paper = 'manylabs'
  write.csv(dplyr::select(tab, experiment, tau2=estimate, ci.lb, ci.ub, paper), paste0('./results/manylabs_vc_exclude_', mm,'.csv'), row.names=F)
}

# ###---Quick thing for Larry
# foo = read.csv("./results/qtest_fixed_manylabs.csv")
# bar = foo %>% 
#   mutate(exact_replication = as.integer(p0 > .05), 
#          approximate_replication = as.integer(p25 > .05)) %>% 
#   select(experiment, k, Q, exact_replication, approximate_replication, mdh_exact = mdh0, mdh_approximate = mdh25)
# 
# bar$power_exact = sapply(experiments, FUN=function(ee){
#       dd = data %>% filter(experiment==ee, is.finite(g))
#       k = nrow(dd)
#       powerRepTest(k=k, v=dd$vg, lambda=k-1, lambda0=0)
#     }
# )
# 
# bar$power_approximate = sapply(experiments, FUN=function(ee){
#   dd = data %>% filter(experiment==ee, is.finite(g))
#   k = nrow(dd)
#   powerRepTest(k=k, v=dd$vg, lambda=k-1, lambda0=(k-1)/4)
# }
# )
# 
# write.csv(bar, "../../paper_psych-methods/manylabs.csv")

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
