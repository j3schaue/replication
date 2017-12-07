###------------------------------------------------------------###
###------------------------------------------------------------###
###------------------------------------------------------------###
### ANALYSES OF RRR:WAGENMAKERS REPLICATION ANALYSES
###------------------------------------------------------------###
###------------------------------------------------------------###
###------------------------------------------------------------###

# set the working directory to src so we can use relative paths
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(dplyr); library(metafor)
source("../package/replicationTest.R")
source("../package/mdh.R")
source("./misc.R")
paper = 'wagenmakers'
tes = 'd'; vr = 'vd'

###------------------------------------------------------------###
###------------------------------------------------------------###
### Wagenmakers Comparison
###------------------------------------------------------------###
###------------------------------------------------------------###
# load the data
data = read.csv("../data/rrr_wagenmakers.csv") %>%
  filter(is.finite(d)) %>%   # weed out infinite estimates
  filter(!is.na(d))

experiments = unique(data$experiment)

df_inc = do.call(rbind, 
                 lapply(experiments, FUN=function(exp){
                   dd = dplyr::filter(data, experiment==exp)
                   ff = rma.uni(dd[[tes]], dd[[vr]], method='FE')
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
                   ff = rma.uni(dd[[tes]], dd[[vr]], method='FE')
                   dd$student = rstudent(ff)$z
                   dd$standard = rstandard(ff)$z
                   dd$Qi = leave1out(ff)$Q
                   dd$tbardot = rep(ff$beta, nrow(dd))
                   return(dd)
                 })
)

##---Sample sizes
foo = df_inc %>% group_by(experiment) %>%
  filter(!(site %in% c('original'))) %>%
  summarize(nmin = min(n), 
            nq1 = quantile(n, .25), 
            nmed = quantile(n, .5),
            nmean = mean(n),
            nq3 = quantile(n, .75), 
            nmax = max(n)) 

##---Check outlier metrics
df_inc %>% 
  group_by(experiment) %>%
  filter(abs(standard) == max(abs(standard)) | Qi == min(Qi)) %>%
  dplyr::select(experiment, site, d, vd, n, standard, Qi)
df_exc %>% 
  group_by(experiment) %>%
  filter(abs(standard) == max(abs(standard)) | Qi == min(Qi)) %>%
  dplyr::select(experiment, site, d, vd, n, standard, Qi)


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
### Wagenmakers Heterogenetiy Q-Test
###------------------------------------------------------------###
###------------------------------------------------------------###


###----Include original study-----------------------------------
# Main analysis
fetab_inc = qtest_results(df_inc, ratios=ratios, t=tes, v=vr, paper=paper, verbose=F) %>% 
  left_join(., distinct(dplyr::select(df_inc, experiment, tbardot)))
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
  left_join(., distinct(dplyr::select(df_exc, experiment, tbardot)))
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
### Wagenmakers variance components
###------------------------------------------------------------###
###------------------------------------------------------------###

# Estimate variance components with Paule-Mandel and DerSimonian-Laird estimators
methods = c("PM", "DL")

###----Include Original Study--------------------------------------------------
for(mm in methods){# for each method, compute tau^2 for each set of replicates
  tab = as.data.frame(do.call(rbind, lapply(experiments, FUN=function(ee){
    dd = filter(df_inc, experiment==ee)
    ff = confint(rma.uni(yi=dd[[tes]], vi=dd[[vr]], method=mm))$random[1,] # get the estimate and the CI
  })))
  tab$experiment = experiments
  tab$vbar = fetab_inc$vbar
  tab$paper = paper
  write.csv(dplyr::select(tab, experiment, tau2=estimate, ci.lb, ci.ub, paper), 
            paste0('./results/', paper, '_vc_include_', mm,'.csv'), 
            row.names=F)
}

###----Exclude Original Study--------------------------------------------------
for(mm in methods){# for each method, compute tau^2 for each set of replicates
  tab = as.data.frame(do.call(rbind, lapply(experiments, FUN=function(ee){
    dd = filter(df_exc, experiment==ee & site!='original') # exclude the original study
    ff = confint(rma.uni(yi=dd[[tes]], vi=dd[[vr]], method=mm))$random[1,] # get the estimate and the CI
  })))
  tab$experiment = experiments
  tab$vbar = fetab_exc$vbar
  tab$paper = paper
  write.csv(dplyr::select(tab, experiment, tau2=estimate, ci.lb, ci.ub, paper),
            paste0('./results/', paper, '_vc_exclude_', mm,'.csv'),
            row.names=F)
}


###------------------------------------------------------------###
###------------------------------------------------------------###
### Wagenmakers Forest Plots
###------------------------------------------------------------###
###------------------------------------------------------------###
plots = list()
for(expt in experiments){
  tmpdf = filter(df_inc, experiment==expt) %>% 
    arrange(as.integer(site != 'original'))
  plots[[expt]] = forest(tmpdf[[tes]], tmpdf[[vr]], slab=tmpdf$site)
  dev.copy(pdf, paste0("./results/plots/", paper, "_", expt, ".pdf"))
  dev.off()
}
