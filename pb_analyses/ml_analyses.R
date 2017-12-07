###------------------------------------------------------------###
###------------------------------------------------------------###
###------------------------------------------------------------###
### ANALYSES OF MANYLABS REPLICATIONS
###------------------------------------------------------------###
###------------------------------------------------------------###
###------------------------------------------------------------###

# set the working directory to src so we can use relative paths
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(dplyr); library(metafor)
source("../package/replicationTest.R")
source("../package/mdh.R")
source("./misc.R")
paper = 'manylabs'
tes = 't'; vr = 'v'

###------------------------------------------------------------###
###------------------------------------------------------------###
### Many Labs Comparison
###------------------------------------------------------------###
###------------------------------------------------------------###
# load the data
data = read.csv("../data/manylabs_comp.csv") %>%
  filter(is.finite(t)) %>%   # weed out infinite estimates
  filter(!is.na(t))

experiments = unique(data$experiment)

df_inc = do.call(rbind, 
             lapply(experiments, FUN=function(exp){
                dd = dplyr::filter(data, experiment==exp)
                ff = rma.uni(dd$t, dd$v, method='FE')
                dd$student = rstudent(ff)$z
                dd$standard = rstandard(ff)$z
                dd$Qi = leave1out(ff)$Q
                dd$tbardot = rep(ff$beta, nrow(dd))
                return(dd)
              })
)


df_exc = do.call(rbind, 
                 lapply(experiments, FUN=function(exp){
                   dd = dplyr::filter(data, experiment==exp & site!='original')
                   ff = rma.uni(dd$t, dd$v, method='FE')
                   dd$student = rstudent(ff)$z
                   dd$standard = rstandard(ff)$z
                   dd$Qi = leave1out(ff)$Q
                   dd$tbardot = rep(ff$beta, nrow(dd))
                   return(dd)
                 })
)

##---Sample sizes
foo = df_inc %>% group_by(experiment) %>%
  filter(!(site %in% c('original', 'mturk', 'pi'))) %>%
  summarize(nmin = min(n), 
            nq1 = quantile(n, .25), 
            nmed = quantile(n, .5),
            nmean = mean(n),
            nq3 = quantile(n, .75), 
            nmax = max(n)) 
write.csv(foo, "mlsamplesizes.csv", row.names=F)

##---Check outlier metrics
df_inc %>% 
  group_by(experiment) %>%
  filter(abs(standard) == max(abs(standard)) | Qi == min(Qi)) %>%
  select(experiment, site, t, v, n, standard, Qi)
df_exc %>% 
  group_by(experiment) %>%
  filter(abs(standard) == max(abs(standard)) | Qi == min(Qi)) %>%
  select(experiment, site, t, v, n, standard, Qi)


#---Comparison Framework--------------------------------
## Take the k-1 replicates and create one synthetic estimator
## using either fixed or random effects meta-analysis.
## Compare the initial finding to this synthetic estimator

# Use a fixed and random effects meta-analysis to combine replications
methods = c('FE', 'DL')

# Set the null hypotheses lambda0 = 0, (k-1)/4, (k-1)/3, 2(k-1)/3
ratios = c(0, 1/4, 1/3, 2/3)

lc = runComparisonAnalyses(data=df_inc, t=tes, v=vr, ratios=ratios, paper=paper, 
                           methods=methods)

for(i in seq(lc)){
  write.csv(lc[[i]],
            paste0("./results/comparison_", paper, "_", names(lc)[i], ".csv"), 
          row.names=F)
}

###------------------------------------------------------------###
###------------------------------------------------------------###
### Many Labs Heterogenetiy Q-Test
###------------------------------------------------------------###
###------------------------------------------------------------###


###----Include original study-----------------------------------
# Main analysis
fetab_inc = qtest_results(df_inc, ratios=ratios, t=tes, v=vr, paper=paper, verbose=F) %>% 
  left_join(., distinct(select(df_inc, experiment, tbardot)))
fetab_inc
write.csv(fetab_inc, 
          paste0("./results/qtest_fixed_", paper, "_include.csv"), row.names=F)

# drop largest outlier
tbardot = do.call(rbind, lapply(experiments, FUN=function(expt){
  dd = df_inc %>% 
    filter(experiment == expt) %>% 
    filter(abs(standard) != max(abs(standard)))
  return(data.frame(experiment = expt, 
                    tbardot = rma.uni(dd[[tes]], dd[[vr]], method='FE')$beta[1,1]))
}))

fetab_inc_ol = qtest_results(df_inc, ratios=ratios, t=tes, v=vr, paper=paper,
                             exclude="abs(standard)!=max(abs(standard))") %>%
  left_join(tbardot)

fetab_inc_ol
write.csv(fetab_inc_ol, 
          paste0("./results/qtest_fixed_", paper, "_include_outlier.csv"), 
          row.names=F)

###----Exclude original study-----------------------------------------------
# Main analysis
fetab_exc = qtest_results(df_exc, ratios=ratios, t=tes, v=vr, paper=paper) %>%
             left_join(., distinct(select(df_exc, experiment, tbardot)))
fetab_exc
write.csv(fetab_exc, 
          paste0("./results/qtest_fixed_", paper, "_exclude.csv"),
          row.names=F)

# drop largest outlier
tbardot = do.call(rbind, lapply(experiments, FUN=function(expt){
  dd = df_exc %>% 
    filter(experiment == expt) %>% 
    filter(abs(standard) != max(abs(standard)))
  return(data.frame(experiment = expt, 
                    tbardot = rma.uni(dd[[tes]], dd[[vr]], method='FE')$beta[1,1]))
}))

fetab_exc_ol = qtest_results(df_exc, ratios=ratios, t=tes, v=vr, paper=paper, 
                             exclude="abs(standard)!=max(abs(standard))") %>%
  left_join(tbardot)

fetab_exc_ol

write.csv(fetab_exc_ol, 
          paste0("./results/qtest_fixed_", paper, "_exclude_outlier.csv"),
          row.names=F)


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
    dd = filter(df_inc, experiment==ee)
    ff = confint(rma.uni(yi=dd$t, vi=dd$v, method=mm))$random[1,] # get the estimate and the CI
  })))
  tab$experiment = experiments
  tab$vbar = fetab_inc$vbar
  tab$paper = 'manylabs'
  write.csv(dplyr::select(tab, experiment, tau2=estimate, ci.lb, ci.ub, paper), 
            paste0('./results/manylabs_vc_include_', mm,'.csv'), 
            row.names=F)
}

###----Exclude Original Study--------------------------------------------------
for(mm in methods){# for each method, compute tau^2 for each set of replicates
  tab = as.data.frame(do.call(rbind, lapply(experiments, FUN=function(ee){
    dd = filter(df_exc, experiment==ee & site!='original') # exclude the original study
    ff = confint(rma.uni(yi=dd$t, vi=dd$v, method=mm))$random[1,] # get the estimate and the CI
  })))
  tab$experiment = experiments
  tab$vbar = fetab_exc$vbar
  tab$paper = 'manylabs'
  write.csv(dplyr::select(tab, experiment, tau2=estimate, ci.lb, ci.ub, paper),
            paste0('./results/', paper, '_vc_exclude_', mm,'.csv'),
            row.names=F)
}


###------------------------------------------------------------###
###------------------------------------------------------------###
### Many Labs Forest Plots
###------------------------------------------------------------###
###------------------------------------------------------------###
plots = list()
for(expt in experiments){
  tmpdf = filter(df_inc, experiment==expt) %>% 
    arrange(as.integer(site != 'original'))
  plots[[expt]] = forest(tmpdf$t, tmpdf$v, slab=tmpdf$site)
  dev.copy(pdf, paste0("./results/plots/", paper, "_", expt, ".pdf"))
  dev.off()
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


# fetab_inc = qtest_results()
# 
# ## Set parameters for analysis
# experiments = unique(df$experiment)
# ks = sapply(experiments, # no. of trials per experiment
#             FUN=function(ee) count(filter(df, experiment==ee))$n)
# tau0s = c(0, 1/4, 1/3, 2/3) # plausible ratios for tau0
# lambda0s = (ks-1)*tau0s  # convert to lambda0
# vbars = sapply(experiments, FUN=function(ee) mean(dplyr::filter(df, experiment==ee)$v)) # avg sampling variances
# fe = lapply(tau0s, FUN=function(tau0) # loop through null hypotheses lambda0 = 0, (k-1)/4, (k-1)/3, 2(k-1)/3
#   setNames(data.frame( # store results as a data frame
#     matrix(unlist(lapply(seq_along(experiments), FUN=function(i)
#       
#       # get the results of the Q-test using all of the studies (rather than aggregating replicates)
#       combineResults(t=filter(df, experiment==experiments[i])$t,
#                      v=filter(df, experiment==experiments[i])$v,
#                      lambda0=(ks[i]-1)*tau0, 
#                      maxratio=20)
#     )
#     ), ncol=5, byrow = T)
#   ), c("k", "Q", paste0("calpha", round(tau0*100, 0)), # name the columns
#        paste0("p", round(tau0*100, 0)), paste0("mdh", round(tau0*100, 0)))
#   )
# )
# 
# # Join for all null hypotheses
# fetab = Reduce(left_join, fe)
# fetab$experiment = experiments
# fetab$vbar = vbars
# fetab$paper = 'manylabs'
# 
# # reorder df and write to file
# fetab = fetab[c("paper", "experiment", "k", "Q", "calpha0",  "p0", "mdh0", "calpha25", "p25",  "mdh25", "calpha33", "p33",  "mdh33", "calpha67", "p67", 
#                 "mdh67", "vbar")] %>%
#   left_join(., distinct(dplyr::select(df, experiment, replicated)))