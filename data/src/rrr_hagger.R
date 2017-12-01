setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")

#Loading Library
library(dplyr)

#Cleaning Data
#RTV
RTV <- read.csv("./raw/rrr_hagger/Data Files and Analysis Files for Meta-Analysis/RTV_incl.csv")
RTV <- select(RTV, c(Study.name, Ego.Depletion.Sample.size, Control.Sample.size, Std.diff.in.means, Std.Err, Hedges.s.g, Std.Err.1))
RTV <- rename(RTV, site = Study.name, n = Ego.Depletion.Sample.size, m = Control.Sample.size, d = Std.diff.in.means, sed = Std.Err, g = Hedges.s.g, seg = Std.Err.1)
RTV["experiment"] <- "RTV"
RTV["es"] <- "d"
RTV["replicated"] <- "0"

#RT
RT <- read.csv("./raw/rrr_hagger/Data Files and Analysis Files for Meta-Analysis/RT_incl.csv")
RT <- select(RT, c(Study.name, Ego.Depletion.Sample.size, Control.Sample.size, Std.diff.in.means, Std.Err, Hedges.s.g, Std.Err.1))
RT <- rename(RT, site = Study.name, n = Ego.Depletion.Sample.size, m = Control.Sample.size, d = Std.diff.in.means, sed = Std.Err, g = Hedges.s.g, seg = Std.Err.1)
RT["experiment"] <- "RT"
RT["es"] <- "d"
RT["replicated"] <- "1"

#Original(Sripada) result
no = 23
mo = 24
j = 1-(3/(4*(no+mo-2)-1))
RTVoriginal <- data.frame(site = "original", n = no, m = mo, d = 0.68, sed = sqrt((no+mo)/(no*mo)+(.68^2)/(2*(no+mo))), g = j*0.68, seg = sqrt((j^2)*(no+mo)/(no*mo)+(.68^2)/(2*(no+mo))), experiment = "RTV", es = "d", replicated = "0")
RToriginal <- data.frame(site = "original", n = no, m = mo, d = .29, sed = sqrt((no+mo)/(no*mo)+(.29^2)/(2*(no+mo))), g = j*0.29, seg = sqrt((j^2)*(no+mo)/(no*mo)+(.29^2)/(2*(no+mo))), experiment = "RT", es = "d", replicated = "1")

#Adding variances for d and g
df <- rbind(RTV, RTVoriginal, RT, RToriginal)
df <- mutate(.data = df,
             vd = sed^2, 
             vg = seg^2) 

write.csv(df,  "rrr_hagger.csv", row.names=F)
