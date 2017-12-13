# set the working directory to src so we can use relative paths
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path), "/results"))
library(dplyr)
source("../misc.R")


###----------------------------------------------------------###
###----------------------------------------------------------###
###----------------------------------------------------------###
### COMPARISON ANALYSES
### Analyses for k=2 comparisons of replication
### Read in analyses on the scale of Cohen's d and
### other scales to check robustness.
###----------------------------------------------------------###
###----------------------------------------------------------###
###----------------------------------------------------------###

#------Check if conversion to Cohen's d affects decisions
# alogna et al
al = read.csv("./comparison_alogna_FE.csv") %>% mutate(es='rd')
ald = read.csv("./comparison_alogna_FE_d.csv") %>% mutate(es='d')
# many labs
ml = read.csv("./comparison_manylabs_FE.csv") %>% 
  filter(experiment %in% c("Allowedforbidden", "Gainloss", "IAT", "Scales")) %>%
  mutate(es = c('logor', 'logor', 'z', 'rd'))
mld = read.csv("./comparison_manylabs_FE_d.csv") %>% 
  filter(experiment %in% c("Allowedforbidden", "Gainloss", "IAT", "Scales")) %>%
  mutate(es='d')
# RPP
rpp = read.csv("./comparison_rpp_FE.csv") %>% mutate(es='z')
rppd = read.csv("./comparison_rpp_FE_d.csv") %>% mutate(es='d')
# RPE
rpe = read.csv("./comparison_rpe_FE.csv") %>% mutate(es='z')
rped = read.csv("./comparison_rpe_FE_d.csv") %>% mutate(es='d')

tocheck = list(alogna=list(al, ald), manylabs=list(ml, mld), rpp=list(rpp, rppd), rpe=list(rpe, rped))
for(i in seq(tocheck)){
  print(names(tocheck)[i])
  print(which(as.integer(tocheck[[i]][[1]]$p0 < .05) != as.integer(tocheck[[i]][[2]]$p0 < .05)))
  print(which(as.integer(tocheck[[i]][[1]]$p25 < .05) != as.integer(tocheck[[i]][[2]]$p25 < .05)))
  print(which(as.integer(tocheck[[i]][[1]]$p33 < .05) != as.integer(tocheck[[i]][[2]]$p33 < .05)))
  print(which(as.integer(tocheck[[i]][[1]]$p67 < .05) != as.integer(tocheck[[i]][[2]]$p67 < .05)))
  
  print(summary(tocheck[[i]][[1]]$p0 - tocheck[[i]][[2]]$p0))
  print(summary(tocheck[[i]][[1]]$p25 - tocheck[[i]][[2]]$p25))
  print(summary(tocheck[[i]][[1]]$p33 - tocheck[[i]][[2]]$p33))
  print(summary(tocheck[[i]][[1]]$p67 - tocheck[[i]][[2]]$p67))
}

# For Alogna et al, there are no differences between conclusions, and p-values differ by less than .02, and on average by .015. 
# Notably for Alogna, all of the (d) p-values are larger than the risk difference p-values
# For Many Labs, there are no differences in conclusions. The largest difference in p-values is 0.03. In fact, the only different
# p-value is for the IAT experiment which was initially on the scale of the z transform. All p-values are about 0.03 larger
# for the z scale. 
# For RPP, there is one experiment that comes to a different conclusion for p25 and p33: E Nurmsoo. It's a small study, and the
# p-values tend to be larger for the z scale.
# There are some larger differences in p-values, but this is largely for smaller studies.
# For RPE, there are no qualitative differences. p-values are all within .02 of each other. 

###------------Cohen's d scale analysis
dtoload = grep("(alogna|manylabs|rpe|rpp)_FE.csv", grep("comparison_.*_FE", list.files(), value=T), value=T, invert=T)
ddfs = lapply(dtoload, FUN=function(ff){
  print(ff)
  foo = read.csv(ff) 
  foo = foo %>%
    dplyr::select(paper, experiment, k, Q, calpha0, p0, mdh0, calpha25, 
                  p25, mdh25, calpha33, p33, mdh33, calpha67, p67, 
                  mdh67, replicated, t1, t2, v1, v2) %>%
    filter(!is.na(Q)) %>%
    mutate(mdh0scale = sqrt(mdh0 * (v1 + v2)), # MDH on the scale of absolute differences
           mdh25scale = sqrt(mdh25 * (v1 + v2)),
           mdh33scale = sqrt(mdh33 * (v1 + v2)),
           mdh67scale = sqrt(mdh67 * (v1 + v2)))
  return(foo)
  
})
dft = do.call(rbind, ddfs)
write.csv(dft, './aggregate/full_comparison_d.csv', row.names=F)

##----aggregate by paper
# group all registered replication reports
rrrs = c('alogna', 'cheung', 'eerland', 'hagger', 'rand', 'wagenmakers')
daggrrr = dft %>% filter(paper %in% rrrs) %>%
  summarize(nstudies = n(), 
            nonreplication = sum(replicated == 0),
            nonrep0 = sum(p0 < .05), # original nonreplication determinations
            nonrep25 = sum(p25 < .05), # nonreplications according to Q test
            nonrep33 = sum(p33 < .05),
            nonrep67 = sum(p67 < .05), 
            mdh0scale = mean(mdh0scale), # MDH on the scale of absolute differences
            mdh25scale = mean(mdh25scale),
            mdh33scale = mean(mdh33scale),
            mdh67scale = mean(mdh67scale)) %>%
  cbind(data.frame(paper = 'rrr'), .)

dagg1 = dft %>% filter(!(paper %in% rrrs)) %>% group_by(paper) %>%
  summarize(nstudies = n(), 
            nonreplication = sum(replicated == 0),
            nonrep0 = sum(p0 < .05), # original nonreplication determinations
            nonrep25 = sum(p25 < .05), # nonreplications according to Q test
            nonrep33 = sum(p33 < .05),
            nonrep67 = sum(p67 < .05), 
            mdh0scale = mean(mdh0scale), # MDH on the scale of absolute differences
            mdh25scale = mean(mdh25scale),
            mdh33scale = mean(mdh33scale),
            mdh67scale = mean(mdh67scale))

dagg = rbind(daggrrr, dagg1)
dagg
write.csv(dagg, "./aggregate/table2d.csv", row.names=F)

##----get most powerful studies
sen = dft %>% filter(mdh0scale <= .3)
sen
write.csv(sen, "./aggregate/send.csv", row.names=F)

# Only 3 studies are powered to detect anything on the oder of an absolute difference smaller than d = .3
# Two are Many Labs studies where the initial variance is quite small (v1 = .0025 for both Allowed/Forbidden and Reciprocity)
# One is an RPP where v1 = .007 and v2 = .001. 
# for what it's worth that's about 1000 people per arm!


###-----Comparison analysis with raw effect sizes-------------------
# Mostly checking results above
toload = grep("(alogna|manylabs|rpe|rpp)_FE_d.csv", grep("comparison_.*_FE", list.files(), value=T), value=T, invert=T)
dfs = lapply(toload, FUN=function(ff){
  print(ff)
  foo = read.csv(ff) 
  foo = foo %>%
    dplyr::select(paper, experiment, k, Q, calpha0, p0, mdh0, calpha25, 
                  p25, mdh25, calpha33, p33, mdh33, calpha67, p67, 
                  mdh67, replicated, t1, t2, v1, v2) %>%
    filter(!is.na(Q)) %>%
    mutate(mdh0scale = sqrt(mdh0 * (v1 + v2)), # MDH on the scale of absolute differences
           mdh25scale = sqrt(mdh25 * (v1 + v2)),
           mdh33scale = sqrt(mdh33 * (v1 + v2)),
           mdh67scale = sqrt(mdh67 * (v1 + v2)))
  return(foo)
  
})
ft = do.call(rbind, dfs)
write.csv(ft, './aggregate/full_comparison.csv', row.names=F)

##----aggregate by paper
# group all registered replication reports
aggrrr = ft %>% filter(paper %in% rrrs) %>%
  summarize(nstudies = n(), 
            nonreplication = sum(replicated == 0),
            nonrep0 = sum(p0 < .05), # original nonreplication determinations
            nonrep25 = sum(p25 < .05), # nonreplications according to Q test
            nonrep33 = sum(p33 < .05),
            nonrep67 = sum(p67 < .05), 
            mdh0scale = mean(mdh0scale), # MDH on the scale of absolute differences
            mdh25scale = mean(mdh25scale),
            mdh33scale = mean(mdh33scale),
            mdh67scale = mean(mdh67scale)) %>%
  cbind(data.frame(paper = 'rrr'), .)

agg1 = ft %>% filter(!(paper %in% rrrs)) %>% group_by(paper) %>%
  summarize(nstudies = n(), 
            nonreplication = sum(replicated == 0),
            nonrep0 = sum(p0 < .05), # original nonreplication determinations
            nonrep25 = sum(p25 < .05), # nonreplications according to Q test
            nonrep33 = sum(p33 < .05),
            nonrep67 = sum(p67 < .05), 
            mdh0scale = mean(mdh0scale), # MDH on the scale of absolute differences
            mdh25scale = mean(mdh25scale),
            mdh33scale = mean(mdh33scale),
            mdh67scale = mean(mdh67scale))

agg = rbind(aggrrr, agg1)
agg
dagg # confirmed it's only the one RPP experiment that has a few different conclusions. 

write.csv(agg, './aggregate/table2.csv', row.names=F)


###----------------------------------------------------------------------###
###----------------------------------------------------------------------###
### Q-test Results
###----------------------------------------------------------------------###
###----------------------------------------------------------------------###

###-------------INCLUDE INITIAL FINDING---------------------------------###

###--------Cohen's d Scale 
dtoloadq = grep("outlier", #drop outlier files
  grep("(alogna|manylabs)_include.csv", # only include the _d files for algona/manylabs 
     grep("qtest_fixed_.*_include", list.files(), value=T), 
     value=T, invert=T),
  value=T, invert=T)
dqs = do.call(rbind, lapply(dtoloadq, read.csv)) %>%
  mutate(mdh0scale = sqrt(2 * vbar * mdh0 / (k - 1)),
         mdh25scale = sqrt(2 * vbar * mdh25 / (k - 1)),
         mdh33scale = sqrt(2 * vbar * mdh33 / (k - 1)),
         mdh67scale = sqrt(2 * vbar * mdh67 / (k - 1)))
dqs
write.csv(dqs, "./aggregate/qtest_include_full_d.csv", row.names=F)

write.csv(dqs %>% select(paper, experiment, p0, p25, p33, p67, mdh0scale, mdh25scale, mdh33scale, mdh67scale), 
          "./aggregate/table3d.csv")

###--------Natural Scale 
toloadq = grep("outlier", #drop outlier files
               grep("(alogna|manylabs)_include_d.csv", # only include the _d files for algona/manylabs 
                    grep("qtest_fixed_.*_include", list.files(), value=T), 
                    value=T, invert=T),
               value=T, invert=T)
qs = do.call(rbind, lapply(toloadq, read.csv)) %>%
  mutate(mdh0scale = sqrt(2 * vbar * mdh0 / (k - 1)),
         mdh25scale = sqrt(2 * vbar * mdh25 / (k - 1)),
         mdh33scale = sqrt(2 * vbar * mdh33 / (k - 1)),
         mdh67scale = sqrt(2 * vbar * mdh67 / (k - 1)))
qs
write.csv(qs, "./aggregate/qtest_include_full.csv", row.names=F)


###----Compare scales
# qs$experiment == dqs$experiment
sum(as.integer(qs$p0 < .05) != as.integer(dqs$p0 < .05))
sum(as.integer(qs$p25 < .05) != as.integer(dqs$p25 < .05))
sum(as.integer(qs$p33 < .05) != as.integer(dqs$p33 < .05))
sum(as.integer(qs$p67 < .05) != as.integer(dqs$p67 < .05))
summary(abs(qs$p0 - dqs$p0))
summary(abs(qs$p25 - dqs$p25))
summary(abs(qs$p33 - dqs$p33))
summary(abs(qs$p67 - dqs$p67))
# p-values are pretty close to the same. Only a few have a difference large enough to care about, but they don't affect conclusions.


###-------------EXCLUDE INITIAL FINDING---------------------------------###

###--------Cohen's d Scale 
dtoloadqexc = grep("outlier", #drop outlier files
                   grep("(alogna|manylabs)_exclude.csv", # only include the _d files for algona/manylabs 
                        grep("qtest_fixed_.*_exclude", list.files(), value=T), 
                        value=T, invert=T),
                   value=T, invert=T)
dqsex = do.call(rbind, lapply(dtoloadqexc, read.csv)) %>%
  mutate(mdh0scale = sqrt(2 * vbar * mdh0 / (k - 1)),
         mdh25scale = sqrt(2 * vbar * mdh25 / (k - 1)),
         mdh33scale = sqrt(2 * vbar * mdh33 / (k - 1)),
         mdh67scale = sqrt(2 * vbar * mdh67 / (k - 1)))
dqsex
write.csv(dqs, "./aggregate/qtest_exclude_full_d.csv", row.names=F)

write.csv(dqs %>% select(paper, experiment, p0, p25, p33, p67, mdh0scale, mdh25scale, mdh33scale, mdh67scale), 
          "./aggregate/table3d_exc.csv")

###--------Natural Scale 
toloadqexc = grep("outlier", #drop outlier files
                  grep("(alogna|manylabs)_exclude_d.csv", # only include the _d files for algona/manylabs 
                       grep("qtest_fixed_.*_exclude", list.files(), value=T), 
                       value=T, invert=T),
                  value=T, invert=T)
qsex = do.call(rbind, lapply(toloadqexc, read.csv)) %>%
  mutate(mdh0scale = sqrt(2 * vbar * mdh0 / (k - 1)),
         mdh25scale = sqrt(2 * vbar * mdh25 / (k - 1)),
         mdh33scale = sqrt(2 * vbar * mdh33 / (k - 1)),
         mdh67scale = sqrt(2 * vbar * mdh67 / (k - 1)))
qsex
write.csv(qsex, "./aggregate/qtest_exclude_full.csv", row.names=F)

###----Compare scales
# qsex$experiment == dqsex$experiment
sum(as.integer(qsex$p0 < .05) != as.integer(dqsex$p0 < .05))
sum(as.integer(qsex$p25 < .05) != as.integer(dqsex$p25 < .05))
sum(as.integer(qsex$p33 < .05) != as.integer(dqsex$p33 < .05))
sum(as.integer(qsex$p67 < .05) != as.integer(dqsex$p67 < .05))
summary(abs(qsex$p0 - dqsex$p0))
summary(abs(qsex$p25 - dqsex$p25))
summary(abs(qsex$p33 - dqsex$p33))
summary(abs(qsex$p67 - dqsex$p67))
# Not sensitive to change in scale.

###----------------------------------------------------------------------###
###----------------------------------------------------------------------###
### Outlier Results
###----------------------------------------------------------------------###
###----------------------------------------------------------------------###

###-------------INCLUDE INITIAL FINDING---------------------------------###

###--------Cohen's d Scale 
dtoloadol = grep("outlier", #drop outlier files
                 grep("(alogna|manylabs)_include_outlier.csv", # only include the _d files for algona/manylabs 
                      grep("qtest_fixed_.*_include", list.files(), value=T), 
                      value=T, invert=T),
                 value=T)
dol = do.call(rbind, lapply(dtoloadol, read.csv)) %>%
  select(paper, experiment, k, outliersite, stdresid, Qi, 
         p0, p25, p33, p67, tbardot) %>%
  left_join(dplyr::select(dqs, paper, experiment, Q)) %>%
  select(paper, experiment, outliersite, stdresid, Q, Qi,
         p0, p25, p33, p67, tbardot)

# Compare outlier analysis to full analysis on Cohen's d scale
# dol$experiment == dqs$experiment
dp0 = which(as.integer(dol$p0 < .05) != as.integer(dqs$p0 < .05))
dp25 = which(as.integer(dol$p25 < .05) != as.integer(dqs$p25 < .05))
dp33 = which(as.integer(dol$p33 < .05) != as.integer(dqs$p33 < .05))
dp67 = which(as.integer(dol$p67 < .05) != as.integer(dqs$p67 < .05))

dol$diff0 = 0; dol$diff0[dp0] = 1
dol$diff25 = 0; dol$diff25[dp25] = 1
dol$diff33 = 0; dol$diff33[dp33] = 1
dol$diff67 = 0; dol$diff67[dp67] = 1
head(dol, 10)
write.csv(dol, "./aggregate/table4d.csv", row.names=F)

# get only the rows where hypothesis tests change based on one outlier
write.csv(dol %>% slice(unique(c(dp0, dp25, dp33, dp67))), 
          "./aggregate/outliers_d_diffs.csv", row.names=F)


###--------Natural Scale 
toloadol = grep("outlier", #drop outlier files
                grep("(alogna|manylabs)_include_outlier_d.csv", # only include the _d files for algona/manylabs 
                     grep("qtest_fixed_.*_include", list.files(), value=T), 
                     value=T, invert=T),
                value=T)
ol = do.call(rbind, lapply(toloadol, read.csv)) %>%
  select(paper, experiment, k, outliersite, stdresid, Qi, 
         p0, p25, p33, p67, tbardot) %>%
  left_join(dplyr::select(qs, paper, experiment, Q)) %>%
  select(paper, experiment, outliersite, stdresid, Q, Qi,
         p0, p25, p33, p67, tbardot)

ol
write.csv(ol, "./aggregate/table4.csv", row.names=F)


###-------Compare scales
# ol$experiment == dol$experiment
which(ol$outliersite != dol$outliersite)
which(as.integer(ol$p0 < .05) != as.integer(dol$p0 < .05))
which(as.integer(ol$p25 < .05) != as.integer(dol$p25 < .05))
which(as.integer(ol$p33 < .05) != as.integer(dol$p33 < .05))
which(as.integer(ol$p67 < .05) != as.integer(dol$p67 < .05))
summary(abs(ol$p0 - dol$p0))
summary(abs(ol$p25 - dol$p25))
summary(abs(ol$p33 - dol$p33))
summary(abs(ol$p67 - dol$p67))
# Not sensitive to change of scale.



###----------------------------------------------------------------------###
###----------------------------------------------------------------------###
### Variance Component Estimates
###----------------------------------------------------------------------###
###----------------------------------------------------------------------###

###-----------------------INCLUDE INITIAL STUDY--------------------------###
# Use Cohen's d scale
###-----Paule-Mandel Estimator
dtolaodvcpm = grep("(alogna|manylabs)_vc_include_PM.csv", grep("include_PM", list.files(), value=T), value=T, invert=T)
dpm = do.call(rbind, lapply(dtolaodvcpm, read.csv)) %>% 
  left_join(dplyr::select(dqs, paper, experiment, k, p0, p25, p33, p67, vbar)) %>%
  select(paper, experiment, k, tau2, cilb=ci.lb, ciub=ci.ub, p0, p25, p33, p67, vbar)

dests = dpm %>%
  filter(p33 > .05) # 'replicating studies' are ones who we fail to reject replication for 1/3
nrow(dests)
summary(dests$tau2/dests$vbar)
summary(dests$ciub/dests$vbar)
summary(sqrt(2*dests$tau2))

write.csv(dests, './aggregate/vcs_include_pm.csv', row.names=F)

###-----DerSimonian & Laird Estimator
dtolaodvcdl = grep("(alogna|manylabs)_vc_include_DL.csv", grep("include_DL", list.files(), value=T), value=T, invert=T)
ddl = do.call(rbind, lapply(dtolaodvcdl, read.csv)) %>% 
  left_join(dplyr::select(dqs, paper, experiment, k, p0, p25, p33, p67, vbar)) %>%
  select(paper, experiment, k, tau2, cilb=ci.lb, ciub=ci.ub, p0, p25, p33, p67, vbar)

destsdl = ddl %>%
  filter(p33 > .05) # 'replicating studies' are ones who we fail to reject replication for 1/3
nrow(destsdl)
summary(destsdl$tau2/destsdl$vbar)
summary(destsdl$ciub/destsdl$vbar)
summary(sqrt(2*destsdl$tau2))

write.csv(destsdl, './aggregate/vcs_include_dl.csv', row.names=F)



# # List of files to load for Q-test of multiple studies
# toload = grep("FE", grep("comparison", list.files(), value=T), value=T)
# dfs = lapply(toload, FUN=function(ff){
#         print(ff)
#         foo = read.csv(ff) 
#         foo = foo %>%
#             dplyr::select(experiment, k, Q, calpha0, p0, mdh0, calpha25, 
#                           p25, mdh25, calpha33, p33, mdh33, calpha67, p67, 
#                           mdh67, replicated, t1, t2, v1, v2, paper) %>%
#             filter(!is.na(Q))
#         if(grepl("_d.csv", ff)){}
#         return(foo)
# 
#   })
# 
# fulltab = do.call(rbind, 
#                   lapply(dfs, FUN=function(df) {df %>% mutate(
#                     t195 = max(abs(t1 + 1.96*sqrt(v1)), abs(t1 - 1.96*sqrt(v1))),
#                     t295 = max(abs(t2 + 1.96*sqrt(v2)), abs(t2 - 1.96*sqrt(v2))),
#                     # Scale MDH on the level of absolute differences
#                     mdh0scale = sqrt(mdh0 * (v1 + v2)), 
#                     mdh25scale = sqrt(mdh25 * (v1 + v2)),
#                     mdh33scale = sqrt(mdh33 * (v1 + v2)), 
#                     mdh67scale = sqrt(mdh67 * (v1 + v2)),
#                     mdh0ratio = (mdh0scale/abs(t1)),
#                     mdh25ratio = (mdh25scale/abs(t1)),
#                     mdh33ratio = (mdh33scale/abs(t1)),
#                     mdh67ratio = (mdh67scale/abs(t1)), 
#                     mdh0ratioub = (mdh0scale/abs(t195)),
#                     mdh25ratioub = (mdh25scale/abs(t195)),
#                     mdh33ratioub = (mdh33scale/abs(t195)),
#                     mdh67ratioub = (mdh67scale/abs(t195))) %>%
#                       dplyr::select(paper, experiment, k, Q, 
#                                     p0, p25, p33, p67, 
#                                     mdh0, mdh25, mdh33, mdh67,
#                                     mdh0scale, mdh25scale, mdh33scale, mdh67scale,
#                                     mdh0ratio, mdh25ratio, mdh33ratio, mdh67ratio,
#                                     mdh0ratioub, mdh25ratioub, mdh33ratioub, mdh67ratioub,
#                                     t1, t195, t2, t295, v1, v2, replicated)
#                   })
# )
# 
# str(fulltab)
# head(fulltab)
# write.csv(fulltab, "./aggregate/full_comparison.csv", row.names=F)
# 
# 
# ##----Summarize results by paper-----------------------------------
# agg = fulltab %>% group_by(paper) %>%
#   summarize(nstudies = n(),
#             nonreplication = sum(replicated==0),
#             nonrep0 = sum(p0 < .05), 
#             nonrep25 = sum(p25 < .05),
#             nonrep33 = sum(p33 < .05),
#             nonrep67 = sum(p67 < .05),
#             mdh0scale = mean(mdh0scale),
#             mdh25scale = mean(mdh25scale),
#             mdh33scale = mean(mdh33scale),
#             mdh67scale = mean(mdh67scale),
#             mdh0ratio = mean(mdh0scale/abs(t1)),
#             mdh25ratio = mean(mdh25scale/abs(t1)),
#             mdh33ratio = mean(mdh33scale/abs(t1)),
#             mdh67ratio = mean(mdh67scale/abs(t1)), 
#             mdh0ratioub = mean(mdh0scale/abs(t195)),
#             mdh25ratioub = mean(mdh25scale/abs(t195)),
#             mdh33ratioub = mean(mdh33scale/abs(t195)),
#             mdh67ratioub = mean(mdh67scale/abs(t195)))
# write.csv(agg, "./aggregate/aggregate_compare.csv", row.names=F)
# 
# ## Table 2
# foo = fulltab
# foo$paper = as.character(fulltab$paper)
# foo[foo$paper %in% c('alogna', 'cheung', 'eerland', 'hagger', 'rand', 'wagenmakers'), ]$paper = "RRR"
# out = foo %>% group_by(paper) %>% 
#   summarize(nstudies = n(),
#             nonreplication=sum(replicated==0), 
#             nonrep0 = sum(p0 < .05), 
#             nonrep25 = sum(p25 < .05),
#             nonrep33 = sum(p33 < .05),
#             nonrep67 = sum(p67 < .05),
#             mdh0scale = mean(mdh0scale),
#             mdh25scale = mean(mdh25scale),
#             mdh33scale = mean(mdh33scale),
#             mdh67scale = mean(mdh67scale),
#             mdh0ratio = mean(mdh0scale/abs(t1)),
#             mdh25ratio = mean(mdh25scale/abs(t1)),
#             mdh33ratio = mean(mdh33scale/abs(t1)),
#             mdh67ratio = mean(mdh67scale/abs(t1)),
#             vbar = mean(v1+ v2),)
# out
# write.csv(out, './aggregate/table2.csv')
# 
# ###---Check possible anomalies for MHD ratio 
# sen = fulltab %>% filter(mdh0ratio < 1) # studies that are powered to detect absolute differences smaller than T1
# 
# sen$powerlb = 1 - pnorm(1.96, mean=abs(sen$t1 - 1.64*sqrt(sen$v1)) / sqrt(sen$v1))
# sen$powerlb2 = 1 - pnorm(1.96, mean=abs(sen$t1 - 1.28*sqrt(sen$v1)) / sqrt(sen$v1))
# sen$powerlb3 = 1 - pnorm(1.96, mean=abs(sen$t1 - .67*sqrt(sen$v1)) / sqrt(sen$v1))
# sen$powermed = 1 - pnorm(1.96, mean=abs(sen$t1)/sqrt(sen$v1))
# 
# ###---Table 3
# out3 = sen 
# unique(out3$paper)
# ml = read.csv("../../data/manylabs_comp.csv")
# rpp = read.csv("../../data/rpp.csv")
# rpe = read.csv("../../data/rpe.csv")
# ppir = read.csv("../../data/ppir.csv")
# out3 = out3 %>%  left_join(do.call(rbind, 
#                             lapply(list(ml, rpp, rpe), FUN=function(x)
#                               distinct(dplyr::select(x, experiment, es))))) #%>%
# out3 = out3 %>% dplyr::select(paper, experiment, 
#          mdh0ratio, mdh25ratio, mdh33ratio, mdh67ratio, powerlb3, powermed, t1,
#          p0, p25, p33, p67)
# 
# out3
# write.csv(out3, './aggregate/table3.csv', row.names=F)
# 
# filter(sen, powerlb < .8) %>% 
#   mutate(rep0 = as.integer(p0 > .05), 
#          N = 1/v1 + 3, 
#          pwrfrac = 1 - pnorm(1.96, mean=abs(.75*t1)/sqrt(v1))) %>% 
#   dplyr::select(paper, experiment, N, powermed, pwrfrac, powerlb, powerlb2, powerlb3, t1, t2, v1, v2, mdh0, replicated, rep0)  
# 
# ##----Plots of sensitivity------------------------------------------
# # library(ggplot2)
# # ggplot(fulltab) + geom_point(aes(abs(t1), mdh0scale, color=paper)) + 
# #   geom_abline(slope=1) +
# #   theme(panel.grid.minor = element_blank()) +
# #   labs(x = "|T1|", 
# #        y = "Scaled MDH (0)")
# # ggsave("./plots/mdh_comparison_0.pdf")
# # 
# # ggplot(fulltab) + geom_point(aes(abs(t1), mdh25scale, color=paper)) + 
# #   geom_abline(slope=1) +
# #   theme(panel.grid.minor = element_blank()) +
# #   labs(x = "|T1|", 
# #        y = "Scaled MDH (1/4)")
# # ggsave("./plots/mdh_comparison_25.pdf")
# # 
# # ggplot(fulltab) + geom_point(aes(abs(t1), mdh33scale, color=paper)) + 
# #   geom_abline(slope=1) +
# #   theme(panel.grid.minor = element_blank()) +
# #   labs(x = "|T1|", 
# #        y = "Scaled MDH (1/3)")
# # ggsave("./plots/mdh_comparison_33.pdf")
# # 
# # ggplot(fulltab) + geom_point(aes(abs(t1), mdh67scale, color=paper)) + 
# #   geom_abline(slope=1) +
# #   theme(panel.grid.minor = element_blank()) +
# #   labs(x = "|T1|", 
# #        y = "Scaled MDH (2/3)")
# # ggsave("./plots/mdh_comparison_67.pdf")
# 
# 
# ###--------------------------------------------------------###
# ### Q-test and variance components
# ###--------------------------------------------------------###
# 
# ###----INCLUDE INITIAL STUDY--------------------------------------
# 
# # lists of results files
# toloadq_inc = grep('outlier', grep("qtest_fixed", grep("include", list.files(), value=T), value=T), value=T, invert=T)
# toloadq_incol = grep('outlier', grep("qtest_fixed", grep("include", list.files(), value=T), value=T), value=T)
# toloadtaudl_inc = grep("_DL", grep("vc_include", list.files(), value=T), value=T)
# toloadtaupm_inc = grep("_PM", grep("vc_include", list.files(), value=T), value=T)
# 
# # Q-test results
# qts_inc = do.call(rbind, lapply(toloadq_inc, read.csv)) %>%
#   dplyr::select(paper, experiment, k, Q, p0, p25, p33, p67, 
#          mdh0, mdh25, mdh33, mdh67, vbar, tbardot, replicated) %>%
#   mutate(mdh0scale = sqrt(mdh0 * vbar),
#          mdh25scale = sqrt(mdh25 * vbar),
#          mdh33scale = sqrt(mdh33 * vbar),
#          mdh67scale = sqrt(mdh67 * vbar),
#          mdh0ratio = sqrt(mdh0 * vbar)/tbardot,
#          mdh25ratio = sqrt(mdh25 * vbar)/tbardot,
#          mdh33ratio = sqrt(mdh33 * vbar)/tbardot,
#          mdh67ratio = sqrt(mdh67 * vbar)/tbardot)
# 
# ###---Table 4
# out4 = qts_inc %>%
#   dplyr::select(paper, experiment, k, Q, p0, p25, p33, p67,
#          mdh0scale, mdh25scale, mdh33scale, mdh67scale)
# 
# rrrs = do.call(rbind, 
#                lapply(grep('rrr', list.files('../../data/'), value=T), FUN=function(x)
#                  read.csv(paste0("../../data/", x)) %>% 
#                    dplyr::select(experiment, es) %>% distinct()))
# rrrs$es = 'd'
# write.csv(out4 %>% 
#             left_join(., distinct(
#               dplyr::select(rbind(rrrs, dplyr::select(ml, experiment, es)), 
#                             experiment, es))), 
#           './aggregate/table4.csv', row.names=F)
# 
# ##-## Good to here.
# qts_incol = do.call(rbind, lapply(toloadq_incol, read.csv)) %>%
#   dplyr::select(paper, experiment, outliersite)
# qts_inc = left_join(qts_inc, qts_incol)
# 
# qts_inc %>% group_by(paper) %>% summarize(nstudies = n(), 
#                                           nonrep0 = sum(p0 < .05),
#                                           nonrep25 = sum(p25 < .05),
#                                           nonrep33 = sum(p33 < .05),
#                                           nonrep67 = sum(p67 < .05),
#                                           nonrep = sum(replicated == 0))
# 
# # Variance component estimates (DerSimonian & Laird, Paule & Mandel)
# dls_inc = do.call(rbind, lapply(toloadtaudl_inc, read.csv)) %>% rename(tauDL=tau2)
# pms_inc = do.call(rbind, lapply(toloadtaupm_inc, read.csv)) %>% rename(tauPM=tau2)
# 
# # Combine the variance components with Q-test results
# vcs_inc = left_join(pms_inc, dls_inc) %>% 
#   select(paper, experiment, tauDL, tauPM, lb=ci.lb, ub=ci.ub) %>%
#   left_join(qts_inc)
# vcs_inc = cbind(vcs_inc,
#                 data.frame(
#                   t(sapply(1:nrow(vcs_inc), 
#                            FUN=function(i) unlist(qCI(vcs_inc$Q[i], vcs_inc$k[i]))))))
# write.csv(vcs_inc, './aggregate/variance_full_include.csv', row.names=F)
# 
# # Get only the variance components for the 'replicating' studies
# reps_inc = vcs_inc %>% filter(p25 > 0.05) %>%
#   select(paper, experiment, k, Q, tauDL, tauPM, lb, ub, vbar, mdh0, mdh25, 
#          p0, p25, lblambda, ublambda, replicated)
# write.csv(reps_inc, './aggregate/variance_replicates_include.csv', row.names=F)
# 
# # averages by paper
# aggs_inc = reps_inc %>% group_by(paper) %>%
#   summarize(tauDL = mean(tauDL), 
#             tauPM = mean(tauPM),
#             lb = mean(lb), 
#             ub = mean(ub), 
#             vbar = mean(vbar), 
#             mdh0 = mean(mdh0), 
#             mdh25 = mean(mdh25), 
#             lblambda = mean(lblambda),
#             ublambda = mean(ublambda))
# write.csv(aggs_inc, "./aggregate/variance_rep-agg_include.csv", row.names=F)
# 
# 
# 
# ###----EXCLUDE INITIAL STUDY--------------------------------------
# 
# # lists of results files
# toloadq = grep("qtest_fixed", grep("exclude", list.files(), value=T), value=T)
# toloadout = grep('outlier', toloadq, value=T)
# toloadq = setdiff(toloadq, toloadout)
# 
# # Q-test results
# qts = do.call(rbind, lapply(toloadq, read.csv)) %>%
#   dplyr::select(paper, experiment, k, Q, 
#                 p0, p25, p33, p67,
#                 mdh0, mdh25, mdh33, mdh67, 
#                 tbardot, vbar, replicated) %>%
#   # compute standard metrics of sensitivity
#   mutate(mdhtau0 = vbar*mdh0, # tau^2 
#          mdhtau25 = vbar*mdh25,
#          mdhtau33 = vbar*mdh33,
#          mdhtau67 = vbar*mdh67, 
#          mdhpwd0 = 2*vbar*mdh0, # average squared pairwise difference
#          mdhpwd25 = 2*vbar*mdh25,
#          mdhpwd33 = 2*vbar*mdh33,
#          mdhpwd67 = 2*vbar*mdh67) # probable errors # max difference # correlation/confounding
# 
# 
# # Summarize results
# 
# 
# 
# 
# 
# 
# qts0 = do.call(rbind, lapply(toloadout, read.csv)) %>%
#   dplyr::select(paper, experiment, k, Q, 
#                 p0, p25, p33, p67,
#                 mdh0, mdh25, mdh33, mdh67, 
#                 tbardot, vbar, replicated)
# 
# combine_res = data.frame(paper = qts$paper,
#            experiment = qts$experiment,
#            Q = qts$Q,
#            Q0 = qts0$Q,
#            p0 = qts$p0, 
#            p0_0 = qts0$p0,
#            p25 = qts$p25, 
#            p25_0 = qts0$p25,
#            p33 = qts$p33, 
#            p33_0 = qts0$p33,
#            p67 = qts$p67, 
#            p67_0 = qts0$p67)
# 
# diffs = sort(unique(c(which(!((combine_res$p0_0 < .05) == (combine_res$p0 < .05))),
#   which(!((combine_res$p25_0 < .05) == (combine_res$p25 < .05))),
#   which(!((combine_res$p33_0 < .05) == (combine_res$p33 < .05))),
#   which(!((combine_res$p67_0 < .05) == (combine_res$p67 < .05))))))
# 
# write.csv(combine_res[diffs,], "./outlier_exclude.csv")
# 
# 
# 
# 
# 
# 
# 
# toloadtaudl = grep("_DL", grep("vc_exclude", list.files(), value=T), value=T)
# toloadtaupm = grep("_PM", grep("vc_exclude", list.files(), value=T), value=T)
# 
# 
# 
# # Variance component estimates (DerSimonian & Laird, Paule & Mandel)
# dls = do.call(rbind, lapply(toloadtaudl, read.csv)) %>% rename(tauDL=tau2)
# pms = do.call(rbind, lapply(toloadtaupm, read.csv)) %>% rename(tauPM=tau2)
# 
# # Combine the variance components with Q-test results
# vcs = left_join(pms, dls) %>% 
#         select(paper, experiment, tauDL, tauPM, lb=ci.lb, ub=ci.ub) %>%
#         left_join(qts)
# vcs = cbind(vcs,
#       data.frame(
#         t(sapply(1:nrow(vcs), FUN=function(i) unlist(qCI(vcs$Q[i], vcs$k[i]))))))
# write.csv(vcs, './aggregate/variance_full_exclude.csv', row.names=F)
# 
# # Get only the variance components for the 'replicating' studies
# reps = vcs %>% filter(p25 > 0.05) %>%
#   select(paper, experiment, k, Q, tauDL, tauPM, lb, ub, vbar, mdh0, mdh25, 
#          p0, p25, lblambda, ublambda, replicated)
# write.csv(reps, './aggregate/variance_replicates_exclude.csv', row.names=F)
# 
# # averages by paper
# aggs = reps %>% group_by(paper) %>%
#   summarize(tauDL = mean(tauDL), 
#             tauPM = mean(tauPM),
#             lb = mean(lb), 
#             ub = mean(ub), 
#             vbar = mean(vbar), 
#             mdh0 = mean(mdh0), 
#             mdh25 = mean(mdh25), 
#             lblambda = mean(lblambda),
#             ublambda = mean(ublambda))
# write.csv(aggs, "./aggregate/variance_rep-agg_exclude.csv", row.names=F)
# 
# 
# 
