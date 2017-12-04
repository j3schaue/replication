# set the working directory to src so we can use relative paths
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path), "/results"))
library(dplyr)
source("../misc.R")

# List of files to load for Q-test of multiple studies
toload = c(
  grep("rrr", grep("FE", grep("comparison", list.files(), value=T), value=T), 
       value=T, invert=T),
            "comparison_rpp.csv", "comparison_rpe.csv")
dfs = lapply(toload[-5], FUN=function(ff){
        print(ff)
        foo = read.csv(ff) %>%
            dplyr::select(experiment, k, Q, calpha0, p0, mdh0, calpha25, 
                          p25, mdh25, calpha33, p33, mdh33, calpha67, p67, 
                          mdh67, replicated, t1, t2, v1, v2, paper) %>%
            filter(!is.na(Q))
        return(foo)

  })



res = do.call(rbind, dfs)
dplyr::filter(res, paper=='manylabs')
sres = res %>% group_by(paper) %>% 
  summarize(total = n(), nonrep0 = sum(p0 < 0.05), nonrep25 = sum(p25 < 0.05), 
            nonreplication = sum(replicated == 0), 
            mdh0 = mean(mdh0),
            mdh25 = mean(mdh25),
            mdh0scale = sqrt(mean(mdh0) * mean(v1 + v2)), 
            mdh25scale = sqrt(mean(mdh25) * mean(v1 + v2)), 
            t1 = mean(t1, na.rm=T),
            t2 = mean(t2, na.rm=T),
            t195 = max(abs(t1 + 1.95*sqrt(v1)), abs(t1 - 1.96*sqrt(v1))),
            t295 = max(abs(t2 + 1.95*sqrt(v2)), abs(t2 - 1.96*sqrt(v2))), 
            v = mean(v1 + v2), 
            v1 = mean(v1), 
            v2 = mean(v2))

write.csv(sres, "./aggregate/agg_comparison.csv", row.names=F)

fulltab = do.call(rbind, 
  lapply(dfs, FUN=function(df) {df %>% mutate(
    mdh0scale = sqrt(mdh0 * (v1 + v2)), 
    mdh25scale = sqrt(mdh25 * (v1 + v2)),
    t195 = max(abs(t1 + 1.95*sqrt(v1)), abs(t1 - 1.96*sqrt(v1))),
    t295 = max(abs(t2 + 1.95*sqrt(v2)), abs(t2 - 1.96*sqrt(v2)))) %>%
    select(paper, experiment, k, Q, p0, p25, mdh0, mdh25, mdh0scale, mdh25scale,
           t1, t195, t2, t295, v1, v2, replicated)
  })
)

write.csv(fulltab, "./aggregate/full_comparison.csv", row.names=F)



###--------------------------------------------------------###
### Q-test and variance components
###--------------------------------------------------------###

###----EXCLUDE INITIAL STUDY--------------------------------------

# lists of results files
toloadq = grep("qtest_fixed", grep("exclude", list.files(), value=T), value=T)
toloadtaudl = grep("_DL", grep("vc_exclude", list.files(), value=T), value=T)
toloadtaupm = grep("_PM", grep("vc_exclude", list.files(), value=T), value=T)

# Q-test results
qts = do.call(rbind, lapply(toloadq, read.csv)) %>%
  select(paper, experiment, k, Q, p0, p25, mdh0, mdh25, vbar, replicated)

# Variance component estimates (DerSimonian & Laird, Paule & Mandel)
dls = do.call(rbind, lapply(toloadtaudl, read.csv)) %>% rename(tauDL=tau2)
pms = do.call(rbind, lapply(toloadtaupm, read.csv)) %>% rename(tauPM=tau2)

# Combine the variance components with Q-test results
vcs = left_join(pms, dls) %>% 
        select(paper, experiment, tauDL, tauPM, lb=ci.lb, ub=ci.ub) %>%
        left_join(qts)
vcs = cbind(vcs,
      data.frame(
        t(sapply(1:nrow(vcs), FUN=function(i) unlist(qCI(vcs$Q[i], vcs$k[i]))))))
write.csv(vcs, './aggregate/variance_full_exclude.csv', row.names=F)

# Get only the variance components for the 'replicating' studies
reps = vcs %>% filter(p25 > 0.05) %>%
  select(paper, experiment, k, Q, tauDL, tauPM, lb, ub, vbar, mdh0, mdh25, 
         p0, p25, lblambda, ublambda, replicated)
write.csv(reps, './aggregate/variance_replicates_exclude.csv', row.names=F)

# averages by paper
aggs = reps %>% group_by(paper) %>%
  summarize(tauDL = mean(tauDL), 
            tauPM = mean(tauPM),
            lb = mean(lb), 
            ub = mean(ub), 
            vbar = mean(vbar), 
            mdh0 = mean(mdh0), 
            mdh25 = mean(mdh25), 
            lblambda = mean(lblambda),
            ublambda = mean(ublambda))
write.csv(aggs, "./aggregate/variance_rep-agg_exclude.csv", row.names=F)



###----INCLUDE INITIAL STUDY--------------------------------------

# lists of results files
toloadq_inc = grep("_rp", grep("qtest_fixed", grep("include", list.files(), value=T), value=T), invert=T, value=T)
toloadtaudl_inc = grep("_DL", grep("vc_include", list.files(), value=T), value=T)
toloadtaupm_inc = grep("_PM", grep("vc_include", list.files(), value=T), value=T)

# Q-test results
qts_inc = do.call(rbind, lapply(toloadq_inc, read.csv)) %>%
  select(paper, experiment, k, Q, p0, p25, mdh0, mdh25, vbar, replicated)

# Variance component estimates (DerSimonian & Laird, Paule & Mandel)
dls_inc = do.call(rbind, lapply(toloadtaudl_inc, read.csv)) %>% rename(tauDL=tau2)
pms_inc = do.call(rbind, lapply(toloadtaupm_inc, read.csv)) %>% rename(tauPM=tau2)

# Combine the variance components with Q-test results
vcs_inc = left_join(pms_inc, dls_inc) %>% 
  select(paper, experiment, tauDL, tauPM, lb=ci.lb, ub=ci.ub) %>%
  left_join(qts_inc)
vcs_inc = cbind(vcs_inc,
            data.frame(
              t(sapply(1:nrow(vcs_inc), 
                       FUN=function(i) unlist(qCI(vcs_inc$Q[i], vcs_inc$k[i]))))))
write.csv(vcs_inc, './aggregate/variance_full_include.csv', row.names=F)

# Get only the variance components for the 'replicating' studies
reps_inc = vcs_inc %>% filter(p25 > 0.05) %>%
  select(paper, experiment, k, Q, tauDL, tauPM, lb, ub, vbar, mdh0, mdh25, 
         p0, p25, lblambda, ublambda, replicated)
write.csv(reps_inc, './aggregate/variance_replicates_include.csv', row.names=F)

# averages by paper
aggs_inc = reps_inc %>% group_by(paper) %>%
  summarize(tauDL = mean(tauDL), 
            tauPM = mean(tauPM),
            lb = mean(lb), 
            ub = mean(ub), 
            vbar = mean(vbar), 
            mdh0 = mean(mdh0), 
            mdh25 = mean(mdh25), 
            lblambda = mean(lblambda),
            ublambda = mean(ublambda))
write.csv(aggs_inc, "./aggregate/variance_rep-agg_include.csv", row.names=F)
