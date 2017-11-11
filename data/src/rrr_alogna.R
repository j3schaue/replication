###---------------------------------------------------------------###
###---------------------------------------------------------------###
### REGISTERED REPLICATION REPORTS (ALOGNA) cleaning
###   raw data csv's modified from Workbooks on OSF
###---------------------------------------------------------------###
###---------------------------------------------------------------###

library(dplyr); library(tidyr); library(ggplot2)
setwd("~/Documents/replication/replication/data/raw_rrr_alogna/") # jakes' path
# setwd("./data/raw_rrr_alogna/") # relative path

rr1 = read.csv("Table1_mod.csv")
rr2 = read.csv("Table2_mod.csv")
names(rr1) = c("authors", "nv", "cv", "fidv", "npv", 
                          "nc", "cc", "fidc", "npc")
names(rr2) = names(rr1)

rr1$authors = as.character(rr1$authors)
rr2$authors = as.character(rr2$authors)

# log OR
rr1$log_orc = log(
              ((rr1$cv/rr1$nv)/((rr1$nv - rr1$cv)/rr1$nv))/
              ((rr1$cc/rr1$nc)/((rr1$nc - rr1$cc)/rr1$nc)))
rr1$vlog_orc = 1/(rr1$cv) + 1/(rr1$nv - rr1$cv) + 
               1/(rr1$cc) + 1/(rr1$nc - rr1$cc)

rr2$log_orc = log(
  ((rr2$cv/rr2$nv)/((rr2$nv - rr2$cv)/rr2$nv))/
    ((rr2$cc/rr2$nc)/((rr2$nc - rr2$cc)/rr2$nc)))
rr2$vlog_orc = 1/(rr2$cv) + 1/(rr2$nv - rr2$cv) + 
  1/(rr2$cc) + 1/(rr2$nc - rr2$cc)

# risk difference
rr1$rd = rr1$cv/rr1$nv - rr1$cc/rr1$nc
rr1$vrd = rr1$cv*(rr1$nv - rr1$cv)/rr1$nv^3 + rr1$cc*(rr1$nc - rr1$cc)/rr1$nc^3

rr2$rd = rr2$cv/rr2$nv - rr2$cc/rr2$nc
rr2$vrd = rr2$cv*(rr2$nv - rr2$cv)/rr2$nv^3 + rr2$cc*(rr2$nc - rr2$cc)/rr2$nc^3



# rr1$exp = 1:nrow(rr1)
# ggplot(rr1) + geom_point(aes(x=rd, y=exp)) + 
#   geom_segment(aes(x=rd + 2*sqrt(vrd) , xend=rd - 2*sqrt(vrd), y=exp, yend=exp))

# Write to file
write.csv(rr1, "../rrr_alogna_rr1_full.csv", row.names=F)
write.csv(rr2, "../rrr_alogna_rr2_full.csv", row.names=F)

write.csv(select(rr1, authors, rd, vrd, logor=log_orc, vlogor=vlog_orc), 
          "../rrr_alogna_rr1.csv", row.names=F)
write.csv(select(rr2, authors, rd, vrd, logor=log_orc, vlogor=vlog_orc),
          "../rrr_alogna_rr2.csv", row.names=F)

