###---------------------------------------------------------------###
###---------------------------------------------------------------###
### REGISTERED REPLICATION REPORTS (ALOGNA) cleaning
###   raw data csv's modified from Workbooks on OSF
### Outcomes are cell counts by group, and ES is risk difference
### There are two main outcomes, correct ID rates and false ID rates.
### We focus only on correct ID rates.
###---------------------------------------------------------------###
###---------------------------------------------------------------###

library(dplyr); library(tidyr); library(ggplot2)

# set the working directory to src so we can use relative paths
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../raw/rrr_alogna")

###---Raw data
rr1 = read.csv("Table1_mod.csv") # RR1 experiments
rr2 = read.csv("Table2_mod.csv") # RR2 experiments

# rename data frames
names(rr1) = c("site", "nv", "cv", "fidv", "npv", 
                          "nc", "cc", "fidc", "npc")
names(rr2) = names(rr1)

# re-structure data
rr1$site = as.character(rr1$site)
rr2$site = as.character(rr2$site)

###---Code effect sizes (or and risk difference---which is what is reported)
# log OR
rr1$log_or = log(
              ((rr1$cv/rr1$nv)/((rr1$nv - rr1$cv)/rr1$nv))/
              ((rr1$cc/rr1$nc)/((rr1$nc - rr1$cc)/rr1$nc)))
rr1$vlog_or = 1/(rr1$cv) + 1/(rr1$nv - rr1$cv) + 
               1/(rr1$cc) + 1/(rr1$nc - rr1$cc)

rr2$log_or = log(
  ((rr2$cv/rr2$nv)/((rr2$nv - rr2$cv)/rr2$nv))/
    ((rr2$cc/rr2$nc)/((rr2$nc - rr2$cc)/rr2$nc)))
rr2$vlog_or = 1/(rr2$cv) + 1/(rr2$nv - rr2$cv) + 
  1/(rr2$cc) + 1/(rr2$nc - rr2$cc)

# risk difference
rr1$rd = rr1$cv/rr1$nv - rr1$cc/rr1$nc
rr1$vrd = rr1$cv*(rr1$nv - rr1$cv)/rr1$nv^3 + rr1$cc*(rr1$nc - rr1$cc)/rr1$nc^3

rr2$rd = rr2$cv/rr2$nv - rr2$cc/rr2$nc
rr2$vrd = rr2$cv*(rr2$nv - rr2$cv)/rr2$nv^3 + rr2$cc*(rr2$nc - rr2$cc)/rr2$nc^3

###---Code experiment before merging
rr1$experiment = "RR1"
rr2$experiment = "RR2"

# merge data
rr = rbind(rr1, rr2)[c("experiment", "site", "nv", "cv", "fidv", "npv", 
                       "nc", "cc", "fidc", "npc",
                       "log_or", "vlog_or", "rd", "vrd")]

# mark effect size
rr$es = 'rd'

# clean up site names
rr$site = sapply(as.character(rr$site), FUN=function(x) 
  tolower(gsub("Shannon McCoy", "Shannon K. McCoy", strsplit(x, ",")[[1]][1])))
rr$site[grepl('mturk', rr$site)] = "mturk"
rr$site[grepl('original', rr$site)] = "original"
rr$replicated = 1 # authors say both studies replicated.

# Get total sample size
rr$n = rr$nc + rr$nv
summary(rr$n)

# Add columns for Cohen's d and variance
rr$d = rr$log_or * sqrt(3)/pi
rr$vd = rr$vlog_or * 3/pi^2

# Write to file
write.csv(rr, "../../rrr_alogna.csv", row.names=F)

rr %>% group_by(experiment) %>% count()

