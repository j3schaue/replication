# set the working directory to src so we can use relative paths
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path), "/results"))
library(dplyr)
source("../misc.R")


###----------------------------------------------------------###
###----------------------------------------------------------###
###----------------------------------------------------------###
### COMPARISON ANALYSES
###----------------------------------------------------------###
###----------------------------------------------------------###
###----------------------------------------------------------###

# List of files to load for Q-test of multiple studies
toload = grep("FE", grep("comparison", list.files(), value=T), value=T)
dfs = lapply(toload, FUN=function(ff){
        print(ff)
        foo = read.csv(ff) 
        foo = foo %>%
            dplyr::select(experiment, k, Q, calpha0, p0, mdh0, calpha25, 
                          p25, mdh25, calpha33, p33, mdh33, calpha67, p67, 
                          mdh67, replicated, t1, t2, v1, v2, paper) %>%
            filter(!is.na(Q))
        return(foo)

  })

fulltab = do.call(rbind, 
                  lapply(dfs, FUN=function(df) {df %>% mutate(
                    t195 = max(abs(t1 + 1.96*sqrt(v1)), abs(t1 - 1.96*sqrt(v1))),
                    t295 = max(abs(t2 + 1.96*sqrt(v2)), abs(t2 - 1.96*sqrt(v2))),
                    # Scale MDH on the level of absolute differences
                    mdh0scale = sqrt(mdh0 * (v1 + v2)), 
                    mdh25scale = sqrt(mdh25 * (v1 + v2)),
                    mdh33scale = sqrt(mdh33 * (v1 + v2)), 
                    mdh67scale = sqrt(mdh67 * (v1 + v2)),
                    mdh0ratio = (mdh0scale/abs(t1)),
                    mdh25ratio = (mdh25scale/abs(t1)),
                    mdh33ratio = (mdh33scale/abs(t1)),
                    mdh67ratio = (mdh67scale/abs(t1)), 
                    mdh0ratioub = (mdh0scale/abs(t195)),
                    mdh25ratioub = (mdh25scale/abs(t195)),
                    mdh33ratioub = (mdh33scale/abs(t195)),
                    mdh67ratioub = (mdh67scale/abs(t195))) %>%
                      dplyr::select(paper, experiment, k, Q, 
                                    p0, p25, p33, p67, 
                                    mdh0, mdh25, mdh33, mdh67,
                                    mdh0scale, mdh25scale, mdh33scale, mdh67scale,
                                    mdh0ratio, mdh25ratio, mdh33ratio, mdh67ratio,
                                    mdh0ratioub, mdh25ratioub, mdh33ratioub, mdh67ratioub,
                                    t1, t195, t2, t295, v1, v2, replicated)
                  })
)

names(fulltab)
write.csv(fulltab, "./aggregate/full_comparison.csv", row.names=F)


##----Summarize results by paper-----------------------------------
agg = fulltab %>% group_by(paper) %>%
  summarize(nstudies = n(),
         nonreplication = sum(replicated==0),
         nonrep0 = sum(p0 < .05), 
         nonrep25 = sum(p25 < .05),
         nonrep33 = sum(p33 < .05),
         nonrep67 = sum(p67 < .05),
         mdh0scale = mean(mdh0scale),
         mdh25scale = mean(mdh25scale),
         mdh33scale = mean(mdh33scale),
         mdh67scale = mean(mdh67scale),
         mdh0ratio = mean(mdh0scale/abs(t1)),
         mdh25ratio = mean(mdh25scale/abs(t1)),
         mdh33ratio = mean(mdh33scale/abs(t1)),
         mdh67ratio = mean(mdh67scale/abs(t1)), 
         mdh0ratioub = mean(mdh0scale/abs(t195)),
         mdh25ratioub = mean(mdh25scale/abs(t195)),
         mdh33ratioub = mean(mdh33scale/abs(t195)),
         mdh67ratioub = mean(mdh67scale/abs(t195)))
write.csv(agg, "./aggregate/aggregate_compare.csv", row.names=F)

##----Plots of sensitivity------------------------------------------
ggplot(fulltab) + geom_point(aes(abs(t1), mdh0scale, color=paper)) + 
  geom_abline(slope=1) +
  theme(panel.grid.minor = element_blank()) +
  labs(x = "|T1|", 
       y = "Scaled MDH (0)")
ggsave("./plots/mdh_comparison_0.pdf")

ggplot(fulltab) + geom_point(aes(abs(t1), mdh25scale, color=paper)) + 
  geom_abline(slope=1) +
  theme(panel.grid.minor = element_blank()) +
  labs(x = "|T1|", 
       y = "Scaled MDH (1/4)")
ggsave("./plots/mdh_comparison_25.pdf")

ggplot(fulltab) + geom_point(aes(abs(t1), mdh33scale, color=paper)) + 
  geom_abline(slope=1) +
  theme(panel.grid.minor = element_blank()) +
  labs(x = "|T1|", 
       y = "Scaled MDH (1/3)")
ggsave("./plots/mdh_comparison_33.pdf")

ggplot(fulltab) + geom_point(aes(abs(t1), mdh67scale, color=paper)) + 
  geom_abline(slope=1) +
  theme(panel.grid.minor = element_blank()) +
  labs(x = "|T1|", 
       y = "Scaled MDH (2/3)")
ggsave("./plots/mdh_comparison_67.pdf")


###--------------------------------------------------------###
### Q-test and variance components
###--------------------------------------------------------###

###----EXCLUDE INITIAL STUDY--------------------------------------

# lists of results files
toloadq = grep("qtest_fixed", grep("exclude", list.files(), value=T), value=T)
toloadout = grep('outlier', toloadq, value=T)
toloadq = setdiff(toloadq, toloadout)

# Q-test results
qts = do.call(rbind, lapply(toloadq, read.csv)) %>%
  dplyr::select(paper, experiment, k, Q, 
                p0, p25, p33, p67,
                mdh0, mdh25, mdh33, mdh67, 
                tbardot, vbar, replicated) %>%
  mutate(mdhtau0 = vbar*mdh0,
         mdhtau25 = vbar*mdh25,
         mdhtau33 = vbar*mdh33,
         mdhtau67 = vbar*mdh67, 
         mdhpwd0 = 2*vbar*mdh0,
         mdhpwd25 = 2*vbar*mdh25,
         mdhpwd33 = 2*vbar*mdh33,
         mdhpwd67 = 2*vbar*mdh67)


# Summarize results






qts0 = do.call(rbind, lapply(toloadout, read.csv)) %>%
  dplyr::select(paper, experiment, k, Q, 
                p0, p25, p33, p67,
                mdh0, mdh25, mdh33, mdh67, 
                tbardot, vbar, replicated)

combine_res = data.frame(paper = qts$paper,
           experiment = qts$experiment,
           Q = qts$Q,
           Q0 = qts0$Q,
           p0 = qts$p0, 
           p0_0 = qts0$p0,
           p25 = qts$p25, 
           p25_0 = qts0$p25,
           p33 = qts$p33, 
           p33_0 = qts0$p33,
           p67 = qts$p67, 
           p67_0 = qts0$p67)

diffs = sort(unique(c(which(!((combine_res$p0_0 < .05) == (combine_res$p0 < .05))),
  which(!((combine_res$p25_0 < .05) == (combine_res$p25 < .05))),
  which(!((combine_res$p33_0 < .05) == (combine_res$p33 < .05))),
  which(!((combine_res$p67_0 < .05) == (combine_res$p67 < .05))))))

write.csv(combine_res[diffs,], "./outlier_exclude.csv")







toloadtaudl = grep("_DL", grep("vc_exclude", list.files(), value=T), value=T)
toloadtaupm = grep("_PM", grep("vc_exclude", list.files(), value=T), value=T)



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

