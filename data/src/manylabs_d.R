#This script reads in the raw data from the ManyLabs project and returns a cleaned dataframe with standardized effect sizes.

#Loading required libraries
library(dplyr)
library(xlsx)

# set the working directory to src so we can use relative paths
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../raw/manylabs/")

######################################################
### Loading and organizing data for 16 experiments ###
######################################################

###---First order approximation of original data
tmp = read.xlsx('manylabs.xlsx', "Table 2")[2:17, 1:3]
names(tmp) = c('experiment', 'd', 'ci')
tmp$lb = sapply(as.character(tmp$ci), FUN=function(x)
                as.numeric(strsplit(x, ", ")[[1]][1]))
tmp$ub = sapply(as.character(tmp$ci), FUN=function(x)
  as.numeric(strsplit(x, ", ")[[1]][2]))
tmp$d[tmp$d=='na'] = NA
tmp$d = as.numeric(as.character(tmp$d))
tmp$vd_est = ((abs(tmp$d - tmp$lb) +  abs(tmp$d - tmp$ub))/(2*1.96))^2
tmp$site = 'original'

#----After reading the original papers, we can't get an original ES for 
#     'Quote attribution', and we can't get an ES variance for 'Allow/Forbid'
#      and we need to correct the ES for 'Gain/loss' and 'Scales
exps = read.csv('manylabs_experiments.csv')
head(exps)

#---Clean up experiment names from original experiments table
tmp$experiment = sapply(tmp$experiment, FUN=function(x) 
  filter(exps, fullname==as.character(x))$experiment)

#---Set up table to write
origs = left_join(tmp, exps)

# Compute variance of effect sizes from original experiments
origs$vd = (origs$nt + origs$nc)/(origs$nt * origs$nc) + 
              origs$d^2/(2*(origs$nt + origs$nc))
origs$vd[is.na(origs$vd)] = 4/origs[is.na(origs$vd),]$n + 
                          origs[is.na(origs$vd),]$d^2/(2*origs[is.na(origs$vd),]$n)

# Compare with estimates from the provided CIs
round(abs(origs$vd_est - origs$vd), 3)

# If we can't verify variances, then we just use the one backed out from 
# the CI the ManyLabs investigators provided
origs$vd[is.na(origs$vd)] = origs$vd_est[is.na(origs$vd)]


# recode gain/loss
data = round(c(152*c(.72, 1-.72), 155*c(.22, 1-.22)), 0)
origs[origs$experiment =="Gainloss", c("d", "vd")] = 
    c(-log(data[1]*data[4]/(data[2]*data[3])) * sqrt(3)/pi,
      (1/data[1] + 1/data[2] + 1/data[3] + 1/data[4]) * 3/pi^2)

# recode scales
# apply continuity correction since there are several 0s in the replicates
data = round(c(64*c(1-.625, .625), 68*c(.162, 1-.162)), 0) + .5
origs[origs$experiment == "Scales", c("d", "vd")] = 
  c(log(data[1]*data[4]/(data[2]*data[3])) * sqrt(3)/pi,
    (1/data[1] + 1/data[2] + 1/data[3] + 1/data[4]) * 3/pi^2)

# recode IAT correlations
n = 243; rr = .42
origs[origs$experiment == "IAT", c("d", "vd")] = 
  c(2*rr/sqrt(1 - rr^2), 
    4 * ((1 - rr^2)^2/(n -1)) / (1 - rr^2)^3)

# Allow/forbidden original
# https://github.com/ManyLabsOpenScience/ManyLabs1/blob/master/Manylabs_OriginalstudiesESCI.R
# apply continuity correction since there are 0s in the replicates
N = 1300 # estimate from the clever Many Labs folks!
Nnot_allow = round(N * .62) + .5# 806
Nallow = round(N * .21, 0) + .5# 273
Nforbid = round(N * .46) + .5# 598
Nnot_forbid = round(N * .39) + .5 # 507
origs[origs$experiment == 'Allowedforbidden', c('d', 'vd')] = # Cohen's d
  c(log(Nnot_allow*Nnot_forbid / (Nallow * Nforbid)) * sqrt(3)/pi,
    (1/Nnot_allow + 1/Nallow + 1/Nnot_forbid + 1/Nforbid) * 3/pi^2)
  

###---fill in NAs for NT, NC
origs[is.na(origs$nt), c('nt', 'nc')] = origs$n[is.na(origs$nt)]/2
origs$es = 'd'


###################################################
#----REPLICATE RESULTS
###################################################

###----Allow/Forbid
af <- read.xlsx("manylabs.xlsx", "Allowed_forbidden")
af<- af[3:(nrow(af) - 1), 1:(ncol(af) - 2)] #Removing summary rows and notes
names(af) = c("site", "allow", "notallow", "forbid", "notforbid", "excluded", "TestStatistics", "ESmd") #Renaming the variables
af["experiment"] <- "Allowedforbidden" #Adding experiment names
af["es"] <- "d" #Adding original effect size measurements
# Convert logOR to d to compare with original experiment
af <- mutate(.data = af,
             t = log(((notallow + .5) * (notforbid + .5))/((allow + .5) * (forbid + .5))) * sqrt(3)/pi, 
             v = (1/(allow + .5) + 1/(notallow + .5) + 1/(forbid + .5) + 1/(notforbid + .5)) * 3/pi^2, 
             es = rep('d', nrow(af)), 
             nt = allow + notallow,
             nc = forbid + notforbid)


###----Anchoring 1
a1 <- read.xlsx("manylabs.xlsx", "Anchoring1")
a1<- a1[3:nrow(a1), 1:(ncol(a1) - 2)]
names(a1) = c("site", "NHigh", "NControl", "NExcluded", "MeanHigh", "MeanControl", "SDHigh", "SDControl", "tev", "tuev", "dfev", "dfuev", "ESuev", "ESmd")
a1["experiment"] <- "Anchoring1"
a1 <- mutate(.data = a1,
             J = 1 - 3/(4*dfev - 1), 
             t = (MeanHigh - MeanControl)/
               sqrt(((NHigh - 1) * SDHigh^2 + 
                       (NControl - 1) * SDControl^2)/
                      (NHigh + NControl - 2)),
             v = (NHigh + NControl) / (NHigh * NControl) + t^2/(2 * (NHigh + NControl)),
             es = rep('d', nrow(a1))) %>%
      rename(nc = NControl, nt = NHigh)


###----Anchoring 2
a2 <- read.xlsx("manylabs.xlsx", "Anchoring2")
a2<- a2[3:nrow(a2), 1:(ncol(a2) - 2)]
names(a2) = c("site", "NHigh", "NControl", "NExcluded", "MeanHigh", "MeanControl", "SDHigh", "SDControl", "tev", "tuev", "dfev", "dfuev", "ESuev", "ESmd")
a2["experiment"] <- "Anchoring2"
a2 <- mutate(.data = a2,
             J = 1 - 3/(4*dfev - 1), 
             t = (MeanHigh - MeanControl)/
               sqrt(((NHigh - 1) * SDHigh^2 + 
                       (NControl - 1) * SDControl^2)/
                      (NHigh + NControl - 2)), 
             v = (NHigh + NControl) / (NHigh * NControl) + t^2/(2 * (NHigh + NControl)),
             es = rep('d', nrow(a2))) %>%
      rename(nc = NControl, nt = NHigh)



###----Anchoring 3
a3 <- read.xlsx("manylabs.xlsx", "Anchoring3")
a3<- a3[3:nrow(a3), 1:(ncol(a3) - 2)]
names(a3) = c("site", "NHigh", "NControl", "NExcluded", "MeanHigh", "MeanControl", "SDHigh", "SDControl", "tev", "tuev", "dfev", "dfuev", "ESuev", "ESmd")
a3["experiment"] <- "Anchoring3"
a3 <- mutate(.data = a3,
             J = 1 - 3/(4*dfev - 1), 
             t = (MeanHigh - MeanControl)/
               sqrt(((NHigh - 1) * SDHigh^2 + 
                       (NControl - 1) * SDControl^2)/
                      (NHigh + NControl - 2)), 
             v = (NHigh + NControl) / (NHigh * NControl) + t^2/(2 * (NHigh + NControl)), 
             es = rep('d', nrow(a3))) %>%
      rename(nc = NControl, nt = NHigh)



###----Anchoring 4
a4 <- read.xlsx("manylabs.xlsx", "Anchoring4")
a4<- a4[3:nrow(a4), 1:(ncol(a4) - 2)]
names(a4) = c("site", "NHigh", "NControl", "NExcluded", "MeanHigh", "MeanControl", "SDHigh", "SDControl", "tev", "tuev", "dfev", "dfuev", "ESuev", "ESmd")
a4["experiment"] <- "Anchoring4"
a4 <- mutate(.data = a4, 
             J = 1 - 3/(4*dfev - 1), 
             t = (MeanHigh - MeanControl)/
               sqrt(((NHigh - 1) * SDHigh^2 + 
                       (NControl - 1) * SDControl^2)/
                      (NHigh + NControl - 2)), 
             v = (NHigh + NControl) / (NHigh * NControl) + t^2/(2 * (NHigh + NControl)), 
             es = rep('d', nrow(a4))) %>%
      rename(nc = NControl, nt = NHigh)



###----Flag Priming
fp <- read.xlsx("manylabs.xlsx", "Flag Priming")
fp<- fp[4:nrow(fp), 1:(ncol(fp) - 2)]
names(fp) = c("site", "NFlag", "NControl", "NExcluded", "MeanFlag", "MeanControl", "SDFlag", "SDControl", "tev", "tuev", "dfev", "dfuev", "ESuev", "ESmd")
fp["experiment"] <- "Flagpriming"
fp <- mutate(.data = fp, 
             J = 1 - 3/(4*dfev - 1), 
             t = (MeanFlag - MeanControl)/
               sqrt(((NFlag - 1) * SDFlag^2 + 
                       (NControl - 1) * SDControl^2)/ 
                      (NFlag + NControl - 2)), 
             v = (NFlag + NControl) / (NFlag * NControl) + t^2/(2 * (NFlag + NControl)), 
             es = rep('d', nrow(fp))) %>%
      rename(nc = NControl, nt = NFlag)



###----Gain/Loss
gl <- read.xlsx("manylabs.xlsx", "Gain_Loss")
gl<- gl[3:nrow(gl), 1:(ncol(gl) - 2)]
names(gl) = c("site", "NGainNorisk", "NLossNorisk", "NGainRisk", "NLossRisk", "NExcluded", "TestStatistics", "ESmd")
gl["experiment"] <- "Gainloss"
gl["es"] <- "d"
gl <- mutate(.data = gl,
             t = log((NGainNorisk*NLossRisk)/(NGainRisk*NLossNorisk)) * sqrt(3)/pi, 
             v = (1/NGainRisk + 1/NGainNorisk + 1/NLossRisk + 1/NLossNorisk) * 3/pi^2,
             es = rep('d', nrow(gl)), 
             nc = NGainRisk + NGainNorisk, 
             nt = NLossRisk + NLossNorisk)


###----Gambler's Fallacy
gf <- read.xlsx("manylabs.xlsx", "Gambler's Fallacy")
gf<- gf[3:nrow(gf), 1:(ncol(gf) - 2)]
names(gf) = c("site", "NThree6", "NControl", "NExcluded", "MeanThree6", "MeanControl", "SDThree6", "SDControl", "tev", "tuev", "dfev", "dfuev", "ESuev", "ESmd")
gf["experiment"] <- "Gamblersfallacy"
gf <- mutate(.data = gf, 
             J = 1 - 3/(4*dfev - 1), 
             t = (MeanThree6 - MeanControl)/sqrt(((NThree6 - 1) * SDThree6^2 +
                                                    (NControl - 1) * SDControl^2)/
                                                   (NThree6 + NControl - 2)), 
             v = (NThree6 + NControl) / (NThree6 * NControl) + t^2/(2 * (NThree6 + NControl)), 
             es = rep('d', nrow(gf))) %>%
      rename(nc = NControl, nt = NThree6)



###----IAT Correlation
iat <- read.xlsx("manylabs.xlsx", "IAT correlation")
iat<- iat[3:nrow(iat), 1:(ncol(iat) - 3)]
names(iat) = c("site", "N" ,"NExcluded", "r")
iat["experiment"] <- "IAT"
iat <- mutate(.data = iat,
              t = 2*r/sqrt(1 - r^2), # Cohen's d (transform)
              v = 4 * ((1 - r^2)^2/(N - 1)) / (1 - r^2)^3, 
              es = rep('d', nrow(iat)), 
              nc = N/2,
              nt = N/2)


###----Imagined Contact
ic <- read.xlsx("manylabs.xlsx", "Imagined Contact")
ic<- ic[3:nrow(ic), 1:(ncol(ic) - 2)]
names(ic) = c("site", "NContact", "NControl", "NExcluded", "MeanContact", "MeanControl", "SDContact", "SDControl", "tev", "tuev", "dfev", "dfuev", "ESuev", "ESmd")
ic["experiment"] <- "Imaginedcontact"
ic <- mutate(.data = ic, 
             J = 1 - 3/(4*dfev - 1), 
             t = (MeanContact - MeanControl)/
               sqrt(((NContact - 1) * SDContact^2 + 
                       (NControl - 1) * SDControl^2)/
                      (NContact + NControl - 2)), 
             v = (NContact + NControl) / (NContact * NControl) + t^2/(2 * (NContact + NControl)), 
             es = rep('d', nrow(ic))) %>%
      rename(nc = NControl, nt = NContact)



###----Math/Gender
mag <- read.xlsx("manylabs.xlsx", "Math_Art Gender")
mag<- mag[3:nrow(mag), 1:(ncol(mag) - 2)]
names(mag) = c("site", "NFemale", "NControl", "NExcluded", "MeanFemale", "MeanControl", "SDFemale", "SDControl", "tev", "tuev", "dfev", "dfuev", "ESuev", "ESmd")
mag$tev <- as.numeric(mag$tev)
mag["experiment"] <- "Mathartgender"
mag <- mutate(.data = mag, 
              J = 1 - 3/(4*dfev - 1), 
              t = (MeanFemale - MeanControl)/
                sqrt(((NFemale - 1) * SDFemale^2 + 
                        (NControl - 1) * SDControl^2)/
                       (NFemale + NControl - 2)), 
              v = (NFemale + NControl) / (NFemale * NControl) + t^2/(2 * (NFemale + NControl)), 
              es = rep('d', nrow(mag))) %>%
      rename(nc = NControl, nt = NFemale)



###----Currency Priming
mp <- read.xlsx("manylabs.xlsx", "Money Priming")
mp<- mp[3:nrow(mp), 1:(ncol(mp) - 2)]
names(mp) = c("site", "NMoney", "NControl", "NExcluded", "MeanMoney", "MeanControl", "SDMoney", "SDControl", "tev", "tuev", "dfev", "dfuev", "ESuev", "ESmd")
mp["experiment"] <- "Moneypriming"
mp <- mutate(.data = mp, 
             J = 1 - 3/(4*dfev - 1), 
             t = (MeanMoney - MeanControl)/
               sqrt(((NMoney - 1) * SDMoney^2 + 
                       (NControl - 1) * SDControl^2)/
                      (NMoney + NControl - 2)), 
             v = (NMoney + NControl) / (NMoney * NControl) + t^2/(2 * (NMoney + NControl)), 
             es = rep('d', nrow(mp))) %>%
      rename(nc = NControl, nt = NMoney)



###----Quote Attribution
qa <- read.xlsx("manylabs.xlsx", "Quote Attribution")
qa<- qa[3:nrow(qa), 1:(ncol(qa) - 2)]
names(qa) = c("site", "NLiked", "NControl", "NExcluded", "MeanLiked", "MeanControl", "SDLiked", "SDControl", "tev", "tuev", "dfev", "dfuev", "ESuev", "ESmd")
qa["experiment"] <- "Quote Attribution"
qa <- mutate(.data = qa, 
             J = 1 - 3/(4*dfev - 1), 
             t = (MeanLiked - MeanControl)/sqrt(((NLiked - 1) * SDLiked^2 + 
                                                   (NControl - 1) * SDControl^2)/
                                                  (NLiked + NControl - 2)), 
             v = (NLiked + NControl) / (NLiked * NControl) + t^2/(2 * (NLiked + NControl)), 
             es = rep('d', nrow(qa))) %>%
      rename(nc = NControl, nt = NLiked)



###----Reciprocity
r <- read.xlsx("manylabs.xlsx", "Reciprocity")
r<- r[3:nrow(r), 1:(ncol(r) - 2)]
names(r) = c("site", "NFirstYes", "NSecondYes", "NFirstNo", "NSecondNo", "NExcluded", "TestStatistics", "ESmd")
r["experiment"] <- "Reciprocity"
# Convert to d to comapre with original experiment
r <- mutate(.data = r,
            t = sqrt(3)/pi * log((NSecondYes*NFirstNo)/(NSecondNo*NFirstYes)), 
            v = 3/pi^2 * (1/NSecondYes + 1/NSecondNo + 1/NFirstYes + 1/NFirstNo),
            es = rep('d', nrow(r)), 
            nt = NSecondYes + NSecondNo,
            nc = NFirstYes + NFirstNo)

###----Scales
s <- read.xlsx("manylabs.xlsx", "Scales")
s<- s[3:nrow(s), 1:(ncol(s) - 2)]
names(s) = c("site", "NLowLess", "NHighLess", "NLowMore", "NHighMore", "NExcluded", "TestStatistics", "ESmd")
s["experiment"] <- "Scales"
# Convert to RD to compare with original experiment
s <- mutate(.data = s,
            t = log(((NHighMore + .5) * (NLowLess + .5))/((NHighLess + .5) * (NLowMore + .5))) * sqrt(3)/pi, 
            v = (1/(NHighMore + .5) + 1/(NLowMore + .5) + 1/(NHighLess + .5) + 1/(NLowLess + .5)) * 3/pi^2,
            es = rep('d', nrow(s)), 
            nt = NHighMore + NHighLess,
            nc = NLowLess + NLowLess)


###----Sunk Costs
sc <- read.xlsx("manylabs.xlsx", "Sunk Costs")
sc<- sc[3:nrow(sc), 1:(ncol(sc) - 2)]
names(sc) = c("site", "NPaid", "NFree", "NExcluded", "MeanPaid", "MeanFree", "SDPaid", "SDFree", "tev", "tuev", "dfev", "dfuev", "ESuev", "ESmd")
sc["experiment"] <- "Sunkcosts"
sc <- mutate(.data = sc, 
             J = 1 - 3/(4*dfev - 1), 
             t = (MeanPaid - MeanFree)/
               sqrt(((NPaid - 1) * SDPaid^2 + 
                       (NFree - 1) * SDFree^2)/
                      (NPaid + NFree - 2)), 
             v = (NPaid + NFree) / (NPaid * NFree) + t^2/(2 * (NPaid + NFree)), 
             es = rep('d', nrow(sc))) %>%
      rename(nt = NPaid, nc = NFree)



###-------------------------------------------------------------------------------###
### Save Clean Data
###-------------------------------------------------------------------------------###

df_list = list(af, a1, a2, a3, a4, fp, gl, gf, iat, ic, mag, mp, qa, r, s, sc)
dfs = rbind(
  do.call(rbind, lapply(df_list, # replicates
                           FUN=function(x) 
                             dplyr::select(x, experiment, site, t, v, es, nt, nc))),
  origs %>% dplyr::select(experiment, site, t=d, v=vd, es, nt, nc)) %>% # original experiments
  left_join(dplyr::select(origs, experiment, replicated)) %>% # get replication status
  arrange(as.character(experiment), as.character(site)) %>% # group by experiment
  mutate(n=nt + nc) # add sample size

# Looks like we calculated an ES for a lab we should have excluded
dfs %>% filter(!is.na(t)) %>% filter(n < 30)
dfs[dfs$experiment=="Mathartgender" & dfs$site=="qccuny2", c('t', 'v')] = c(NA, NA)

# names(df_list) = unique(dfs$experiment)
# saveRDS(df_list, "../../manylabs_raw.RDS")

str(dfs)
unique(dfs$es)

write.csv(dfs, "../../manylabs_d.csv", row.names=F)

