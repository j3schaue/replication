
# set the working directory to src so we can use relative paths
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

dfs = list(
  rpp = read.csv("./results/qtest_fixed_rpp.csv"),
  rpe = read.csv("./results/qtest_fixed_rpe.csv"),
  manylabsf = read.csv("./results/comparison_manylabs_FE.csv"),
  manylabsr = read.csv("./results/comparison_manylabs_DL.csv"),
  alognaf = read.csv("./results/comparison_alogna_FE.csv"),
  alognar = read.csv("./results/comparison_alogna_DL.csv")
)

lambda0s = c(0, 1/4, 1/3, 2/3)
strlambdas = round(100*lambda0s, 0)

prop_replicated = setNames(data.frame(matrix(unlist(lapply(strlambdas, FUN=function(pp)
  sapply(dfs, FUN=function(df) mean(df[[paste0("p", pp)]] > 0.05, na.rm=T))
)), ncol = length(dfs), byrow=T)), names(dfs))

prop_replicated$lambda0 = lambda0s
prop_replicated

agree = setNames(
  data.frame(
    matrix(unlist(
          lapply(strlambdas, FUN=function(pp)
            sapply(dfs, FUN=function(df){
              mean(df$replicated == as.integer(df[[paste0('p', pp)]] > .05), na.rm=T)
          })
)), ncol=length(dfs), byrow=T)), names(dfs))

agree$lambda0 = lambda0s
agree

# original studies say that findings don't replicate, 
# but the Q test is inconclusive
falsepos = setNames(
  data.frame(
    matrix(unlist(
      lapply(strlambdas, FUN=function(pp)
        sapply(dfs, FUN=function(df){
          mean((df$replicated == 0 & df[[paste0('p', pp)]] > .05), na.rm=T)
        })
      )), ncol=length(dfs), byrow=T)), names(dfs))

falsepos$lambda0 = lambda0s
falsepos

# Original studies say that studies replicate, but
# the Q test disagrees
falseneg = setNames(
  data.frame(
    matrix(unlist(
      lapply(strlambdas, FUN=function(pp)
        sapply(dfs, FUN=function(df){
          mean((df$replicated == 1 & df[[paste0('p', pp)]] < .05), na.rm=T)
        })
      )), ncol=length(dfs), byrow=T)), names(dfs))

falseneg$lambda0 = lambda0s
falseneg

for(pp in strlambdas){
  assign(paste0("disagrees", pp),
         lapply(dfs, FUN=function(df){
            df[which(df$replicated == as.integer(df[[paste0("p", pp)]] < 0.05)), ]
         })
  )
}
  
for(i in 1:length(disagrees0)){
  print(ncol(disagrees0[[i]]))
  if(nrow(disagrees0[[i]]) > 0){
    disagrees0[[i]]$paper = names(disagrees0)[i] 
  }
}

###---------------------------------------------###
# Where do we disagree and why?
###---------------------------------------------###
library(metafor); library(tidyr); library(dplyr)

###---RPP
rpp = read.csv("../data/rpp.csv") %>% select(-exp_name, -X)
disagrees0$rpp %>% select(experiment, k, Q, calpha0, p0, replicated)
rpp$site = gsub("([0-9])", "", gsub("_", "", rpp$site))
head(rpp)
rppz = rpp %>% select(experiment, site, z, pvalo, pvalr) %>%
  spread(site, z) %>% rename(z = orig, zrep=rep)
rpp_comp = rpp %>% select(experiment, site, vz) %>%
  spread(site, vz) %>% rename(v = orig, vrep=rep) %>%
  left_join(rppz) %>%
  left_join(disagrees0$rpp %>% select(experiment, Q, calpha0, p0, replicated, cirep, meta), .) %>%
  mutate(rep_pval = as.integer((pvalo <.05 & pvalr < .05) | (pvalo > .05 & pvalr > .05)))
rpp_comp %>% mutate(pd = abs(z - zrep)/z, qtest = as.integer(p0 > .05)) %>% 
  select(experiment, p0, qtest, replicated, cirep, meta, rep_pval, pd)

###---RPP
rpe = read.csv("../data/rpe.csv") %>% select(-site, -es, -n, -r, -ref) %>% 
  left_join(disagrees0$rpe %>% select(experiment, k, Q, calpha0, p0, replicated), .) %>%
  mutate(stat = z/sqrt(vz))
rpe

###---Many Labs
ml = read.csv("../data/manylabs_comparison.csv")
ml_orig = filter(ml, site=='original') %>% select(experiment, t, v, es, replicated)
ml_reps = data.frame(experiment=NULL, trep=NULL, vrep=NULL, es=NULL, replicated=NULL)
for(ee in unique(ml$experiment)){
  dd = filter(ml, site != 'original' & experiment == ee & !is.infinite(t))
  meta = rma.uni(yi=dd$t, vi=dd$v, method='FE')
  toadd = data.frame(experiment=ee, trep=meta$beta, vrep=meta$se^2, es=unique(dd$es), replicated=unique(dd$replicated))
  ml_reps = rbind(ml_reps, toadd)
}
ml_compare = left_join(ml_reps, ml_orig)

disagrees0$manylabsf %>% select(experiment, k, Q, calpha0, p0, replicated) %>%
  left_join(ml_compare) %>%
  mutate(stat = abs(t)/sqrt(v), statrep = abs(trep)/sqrt(vrep))
