#This script reads in the raw data from the ManyLabs project and returns a cleaned dataframe with standardized effect sizes.

#Loading required libraries
library(dplyr)
library(xlsx)

# set the working directory to src so we can use relative paths
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../raw/manylabs/")

######################################################
### Recode original experiment effect sizes        ###
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
tmp$vd_est = ((tmp$ub - tmp$lb)/(2*1.96))^2
tmp$site = 'original'
tmp$experiment = sapply(tmp$experiment, FUN=function(x) 
  filter(exps, fullname==as.character(x))$experiment)

#----After reading the original papers, we can't get an original ES for 
#     'Quote attribution', and we can't get an ES variance for 'Allow/Forbid'
#      and we need to correct the ES for 'Gain/loss' and 'Scales
exps = read.csv('manylabs_experiments.csv')
head(exps)
origs = left_join(tmp, exps)

# Compute variance of effect sizes from original experiments
origs$vd = (origs$nt + origs$nc)/(origs$nt * origs$nc) + 
  origs$d^2/(2*(origs$nt + origs$nc))
origs$vd[is.na(origs$vd)] = 4/origs[is.na(origs$vd),]$n + 
  origs[is.na(origs$vd),]$d^2/(2*origs[is.na(origs$vd),]$n)
origs$vv = 4/origs$n + origs$d^2/(2*origs$n)

# Compare with estimates from the provided CIs
round(abs(origs$vd_est - origs$vd), 3)
round(abs(origs$vd_est - origs$vv), 3)

# If we can't verify variances, then we just use the one backed out from 
# the CI the ManyLabs investigators provided
origs$vd[is.na(origs$vd)] = origs$vd_est[is.na(origs$vd)]


# recode gain/loss
data = round(c(152*c(.72, 1-.72), 155*c(.22, 1-.22)), 0)
origs[origs$experiment=="Gainloss", c("d", "vd")] = 
  c(-log(data[1]*data[4]/(data[2]*data[3])),
    (1/data[1] + 1/data[2] + 1/data[3] + 1/data[4]))

# recode scales
data = round(c(64*c(1-.625, .625), 68*c(.162, 1-.162)), 0)
origs[origs$experiment=="Scales", c("d", "vd")] = 
  c(data[1]/(data[1] + data[2]) - data[3]/(data[3] + data[4]),
    data[1]*data[2]/(data[1] + data[2])^3 + data[3]*data[4]/(data[3] + data[4])^2)

# recode IAT correlations
n = 213; rr = .42
origs[origs$experiment=="IAT", c("d", "vd")] = 
  c(0.5 * log((1 + rr)/(1 - rr)), 
    vz = 1/(n-3))


######################################################
### Loading and organizing data for 16 experiments ###
######################################################

af <- read.xlsx("manylabs.xlsx", "Allowed_forbidden")
af<- af[3:(nrow(af) - 1), 1:(ncol(af) - 2)] #Removing summary rows and notes
names(af) = c("Site", "NAllowYes", "NAllowNo", "NForbidYes", "NForbidNo", "NExcluded", "TestStatistics", "ESmd") #Renaming the variables
af["Experiment"] <- "Allowedforbidden" #Adding experiment names
af["es.measurement"] <- "md" #Adding original effect size measurements

a1 <- read.xlsx("manylabs.xlsx", "Anchoring1")
a1<- a1[3:nrow(a1), 1:(ncol(a1) - 2)]
names(a1) = c("Site", "NHigh", "NControl", "NExcluded", "MeanHigh", "MeanControl", "SDHigh", "SDControl", "tev", "tuev", "dfev", "dfuev", "ESuev", "ESmd")
a1["Experiment"] <- "Anchoring1"
a1["es.measurement"] <- "md"

a2 <- read.xlsx("manylabs.xlsx", "Anchoring2")
a2<- a2[3:nrow(a2), 1:(ncol(a2) - 2)]
names(a2) = c("Site", "NHigh", "NControl", "NExcluded", "MeanHigh", "MeanControl", "SDHigh", "SDControl", "tev", "tuev", "dfev", "dfuev", "ESuev", "ESmd")
a2["Experiment"] <- "Anchoring2"
a2["es.measurement"] <- "md"

a3 <- read.xlsx("manylabs.xlsx", "Anchoring3")
a3<- a3[3:nrow(a3), 1:(ncol(a3) - 2)]
names(a3) = c("Site", "NHigh", "NControl", "NExcluded", "MeanHigh", "MeanControl", "SDHigh", "SDControl", "tev", "tuev", "dfev", "dfuev", "ESuev", "ESmd")
a3["Experiment"] <- "Anchoring3"
a3["es.measurement"] <- "md"

a4 <- read.xlsx("manylabs.xlsx", "Anchoring4")
a4<- a4[3:nrow(a4), 1:(ncol(a4) - 2)]
names(a4) = c("Site", "NHigh", "NControl", "NExcluded", "MeanHigh", "MeanControl", "SDHigh", "SDControl", "tev", "tuev", "dfev", "dfuev", "ESuev", "ESmd")
a4["Experiment"] <- "Anchoring4"
a4["es.measurement"] <- "md"

fp <- read.xlsx("manylabs.xlsx", "Flag Priming")
fp<- fp[4:nrow(fp), 1:(ncol(fp) - 2)]
names(fp) = c("Site", "NFlag", "NControl", "NExcluded", "MeanFlag", "MeanControl", "SDFlag", "SDControl", "tev", "tuev", "dfev", "dfuev", "ESuev", "ESmd")
fp["Experiment"] <- "Flagpriming"
fp["es.measurement"] <- "md"

gl <- read.xlsx("manylabs.xlsx", "Gain_Loss")
gl<- gl[3:nrow(gl), 1:(ncol(gl) - 2)]
names(gl) = c("Site", "NGainNorisk", "NLossNorisk", "NGainRisk", "NLossRisk", "NExcluded", "TestStatistics", "ESmd")
gl["Experiment"] <- "Gainloss"
gl["es.measurement"] <- "md"

gf <- read.xlsx("manylabs.xlsx", "Gambler's Fallacy")
gf<- gf[3:nrow(gf), 1:(ncol(gf) - 2)]
names(gf) = c("Site", "NThree6", "NControl", "NExcluded", "MeanThree6", "MeanControl", "SDThree6", "SDControl", "tev", "tuev", "dfev", "dfuev", "ESuev", "ESmd")
gf["Experiment"] <- "Gamblersfallacy"
gf["es.measurement"] <- "md"

iat <- read.xlsx("manylabs.xlsx", "IAT correlation")
iat<- iat[3:nrow(iat), 1:(ncol(iat) - 3)]
names(iat) = c("Site", "N" ,"NExcluded", "r")
iat["Experiment"] <- "IAT"
iat["es.measurement"] <- "corr"

ic <- read.xlsx("manylabs.xlsx", "Imagined Contact")
ic<- ic[3:nrow(ic), 1:(ncol(ic) - 2)]
names(ic) = c("Site", "NContact", "NControl", "NExcluded", "MeanContact", "MeanControl", "SDContact", "SDControl", "tev", "tuev", "dfev", "dfuev", "ESuev", "ESmd")
ic["Experiment"] <- "Imaginedcontact"
ic["es.measurement"] <- "md"

mag <- read.xlsx("manylabs.xlsx", "Math_Art Gender")
mag<- mag[3:nrow(mag), 1:(ncol(mag) - 2)]
names(mag) = c("Site", "NFemale", "NControl", "NExcluded", "MeanFemale", "MeanControl", "SDFemale", "SDControl", "tev", "tuev", "dfev", "dfuev", "ESuev", "ESmd")
mag$tev <- as.numeric(mag$tev)
mag["Experiment"] <- "Mathartgender"
mag["es.measurement"] <- "md"

mp <- read.xlsx("manylabs.xlsx", "Money Priming")
mp<- mp[3:nrow(mp), 1:(ncol(mp) - 2)]
names(mp) = c("Site", "NMoney", "NControl", "NExcluded", "MeanMoney", "MeanControl", "SDMoney", "SDControl", "tev", "tuev", "dfev", "dfuev", "ESuev", "ESmd")
mp["Experiment"] <- "Moneypriming"
mp["es.measurement"] <- "md"

qa <- read.xlsx("manylabs.xlsx", "Quote Attribution")
qa<- qa[3:nrow(qa), 1:(ncol(qa) - 2)]
names(qa) = c("Site", "NLiked", "NControl", "NExcluded", "MeanLiked", "MeanControl", "SDLiked", "SDControl", "tev", "tuev", "dfev", "dfuev", "ESuev", "ESmd")
qa["Experiment"] <- "Quote Attribution"
qa["es.measurement"] <- "md"

r <- read.xlsx("manylabs.xlsx", "Reciprocity")
r<- r[3:nrow(r), 1:(ncol(r) - 2)]
names(r) = c("Site", "NFirstYes", "NSecondYes", "NFirstNo", "NSecondNo", "NExcluded", "TestStatistics", "ESmd")
r["Experiment"] <- "Reciprocity"
r["es.measurement"] <- "md"

s <- read.xlsx("manylabs.xlsx", "Scales")
s<- s[3:nrow(s), 1:(ncol(s) - 2)]
names(s) = c("Site", "NLowLess", "NHighLess", "NLowMore", "NHighMore", "NExcluded", "TestStatistics", "ESmd")
s["Experiment"] <- "Scales"
s["es.measurement"] <- "md"

sc <- read.xlsx("manylabs.xlsx", "Sunk Costs")
sc<- sc[3:nrow(sc), 1:(ncol(sc) - 2)]
names(sc) = c("Site", "NPaid", "NFree", "NExcluded", "MeanPaid", "MeanFree", "SDPaid", "SDFree", "tev", "tuev", "dfev", "dfuev", "ESuev", "ESmd")
sc["Experiment"] <- "Sunkcosts"
sc["es.measurement"] <- "md"

#Test to see why their reported effect size measurements were inaccurate

#############################################################################
###For experiments with categorical results, generate log odds ratio first###
#############################################################################

af <- mutate(.data = af,
             oddsratio = (NAllowYes*NForbidNo)/(NAllowNo*NForbidYes), #Odds ratio
             log.oddsratio = log(oddsratio), #Log odds ratio
             var.log.oddsratio = 1/NAllowNo + 1/NAllowYes + 1/NForbidYes + 1/NForbidNo, #Variance of log odds ratio
             d = log.oddsratio*(sqrt(3)/pi), #Cohen's D
             g = d*(1 - 3/(4*(NAllowNo + NAllowYes + NForbidYes + NForbidNo - 2) - 1)), #Hedges' G
             vd = var.log.oddsratio*(3/pi^2), #Variance of D
             vg = vd*((1 - 3/(4*(NAllowNo + NAllowYes + NForbidYes + NForbidNo - 2)-1))^2), #Variance of G
             a = (((NAllowYes + NAllowNo)+(NForbidYes + NForbidNo))^2)/((NAllowYes + NAllowNo)*(NForbidYes + NForbidNo)), #Correction factor
             r = d/sqrt(d^2+a),#Correlation
             var.r = ((a^2)*vd)/(((d^2)+a)^3)) #Variance of correlation


gl <- mutate(.data = gl,
             oddsratio = (NGainRisk*NLossNorisk)/(NGainNorisk*NLossRisk),
             log.oddsratio = log(oddsratio), 
             var.log.oddsratio = 1/NGainRisk + 1/NGainNorisk + 1/NLossRisk + 1/NLossNorisk,
             d = log.oddsratio*(sqrt(3)/pi),
             g = d*(1 - 3/(4*(NGainRisk + NGainNorisk + NLossRisk + NLossNorisk - 2) - 1)),
             vd = var.log.oddsratio*(3/pi^2),
             vg = vd*((1 - 3/(4*(NGainRisk + NGainNorisk + NLossRisk + NLossNorisk - 2)-1))^2),
             a = (((NGainRisk + NGainNorisk)+(NLossRisk + NLossNorisk))^2)/((NGainRisk + NGainNorisk)*(NLossRisk + NLossNorisk)),
             r = d/sqrt(d^2+a),
             var.r = ((a^2)*vd)/(((d^2)+a)^3))

r <- mutate(.data = r,
            oddsratio = (NSecondYes*NFirstNo)/(NSecondNo*NFirstYes),
            log.oddsratio = log(oddsratio), 
            var.log.oddsratio = 1/NSecondYes + 1/NSecondNo + 1/NFirstYes + 1/NFirstNo,
            d = log.oddsratio*(sqrt(3)/pi),
            g = d*(1 - 3/(4*(NFirstYes + NSecondYes + NFirstNo + NSecondNo - 2) - 1)),
            vd = var.log.oddsratio*(3/pi^2),
            vg = vd*((1 - 3/(4*(NFirstYes + NSecondYes + NFirstNo + NSecondNo - 2)-1))^2),
            a = (((NSecondYes + NSecondNo)+(NFirstYes + NFirstNo))^2)/((NSecondYes + NSecondNo)*(NFirstYes + NFirstNo)),
            r = d/sqrt(d^2+a),
            var.r = ((a^2)*vd)/(((d^2)+a)^3))

s <- mutate(.data = s,
            oddsratio = (NHighMore*NLowLess)/(NHighLess*NLowMore),
            log.oddsratio = log(oddsratio), 
            var.log.oddsratio = 1/NHighMore + 1/NHighLess + 1/NLowMore + 1/NLowLess,
            d = log.oddsratio*(sqrt(3)/pi),
            g = d*(1 - 3/(4*(NLowLess + NHighLess + NLowMore + NHighMore - 2) - 1)),
            vd = var.log.oddsratio*(3/pi^2),
            vg = vd*((1 - 3/(4*(NLowLess + NHighLess + NLowMore + NHighMore - 2)-1))^2),
            a = (((NHighMore + NHighLess)+(NLowMore + NLowLess))^2)/((NHighMore + NHighLess)*(NLowMore + NLowLess)),
            r = d/sqrt(d^2+a),
            var.r = ((a^2)*vd)/(((d^2)+a)^3))

#####################################################################################
###For experiments with continuous results, generate the mean difference (d) first###
#####################################################################################

a1 <- mutate(.data = a1,
             pooled_var = ((NHigh - 1) * SDHigh^2 + (NControl - 1) * SDControl^2)/(NHigh + NControl - 2), # pooled variance
             J = 1 - 3/(4*dfev - 1), # small sample correction for G
             d = (MeanHigh - MeanControl)/sqrt(pooled_var), #Cohen's D
             g = d * J, #Hedges' G
             vd = (NHigh + NControl) / (NHigh * NControl) + d^2/(2 * (NHigh + NControl)), # Variance of D
             vg = vd * J^2, # Variance of G
             log.oddsratio = d*(pi/sqrt(3)),#Log odds ratio
             var.log.oddsratio = vd*(pi^2)/3, #Variance of log odds ratio
             a = ((NHigh+NControl)^2)/(NHigh*NControl), #Correction factor
             r = d/sqrt(d^2+a), #Correlation
             var.r = ((a^2)*vd)/(((d^2)+a)^3)) #Variance of correlation

a2 <- mutate(.data = a2,
             pooled_var = ((NHigh - 1) * SDHigh^2 + (NControl - 1) * SDControl^2)/(NHigh + NControl - 2), 
             J = 1 - 3/(4*dfev - 1), 
             d = (MeanHigh - MeanControl)/sqrt(pooled_var), 
             g = d * J, 
             vd = (NHigh + NControl) / (NHigh * NControl) + d^2/(2 * (NHigh + NControl)), 
             vg = vd * J^2, 
             log.oddsratio = d*(pi/sqrt(3)),
             var.log.oddsratio = vd*(pi^2)/3,
             a = ((NHigh+NControl)^2)/(NHigh*NControl),
             r = d/sqrt(d^2+a),
             var.r = ((a^2)*vd)/(((d^2)+a)^3)) 

a3 <- mutate(.data = a3,
             pooled_var = ((NHigh - 1) * SDHigh^2 + (NControl - 1) * SDControl^2)/(NHigh + NControl - 2), 
             J = 1 - 3/(4*dfev - 1), 
             d = (MeanHigh - MeanControl)/sqrt(pooled_var), 
             g = d * J, 
             vd = (NHigh + NControl) / (NHigh * NControl) + d^2/(2 * (NHigh + NControl)), 
             vg = vd * J^2, 
             log.oddsratio = d*(pi/sqrt(3)),
             var.log.oddsratio = vd*(pi^2)/3,
             a = ((NHigh+NControl)^2)/(NHigh*NControl),
             r = d/sqrt(d^2+a),
             var.r = ((a^2)*vd)/(((d^2)+a)^3)) 


a4 <- mutate(.data = a4,
             pooled_var = ((NHigh - 1) * SDHigh^2 + (NControl - 1) * SDControl^2)/(NHigh + NControl - 2), 
             J = 1 - 3/(4*dfev - 1), 
             d = (MeanHigh - MeanControl)/sqrt(pooled_var), 
             g = d * J, 
             vd = (NHigh + NControl) / (NHigh * NControl) + d^2/(2 * (NHigh + NControl)), 
             vg = vd * J^2, 
             log.oddsratio = d*(pi/sqrt(3)),
             var.log.oddsratio = vd*(pi^2)/3,
             a = ((NHigh+NControl)^2)/(NHigh*NControl),
             r = d/sqrt(d^2+a),
             var.r = ((a^2)*vd)/(((d^2)+a)^3)) 

fp <- mutate(.data = fp,
             pooled_var = ((NFlag - 1) * SDFlag^2 + (NControl - 1) * SDControl^2)/(NFlag + NControl - 2), 
             J = 1 - 3/(4*dfev - 1), 
             d = (MeanFlag - MeanControl)/sqrt(pooled_var), 
             g = d * J, 
             vd = (NFlag + NControl) / (NFlag * NControl) + d^2/(2 * (NFlag + NControl)), 
             vg = vd * J^2, 
             log.oddsratio = d*(pi/sqrt(3)),
             var.log.oddsratio = vd*(pi^2)/3,
             a = ((NFlag+NControl)^2)/(NFlag*NControl),
             r = d/sqrt(d^2+a),
             var.r = ((a^2)*vd)/(((d^2)+a)^3)) 

gf <- mutate(.data = gf,
             pooled_var = ((NThree6 - 1) * SDThree6^2 + (NControl - 1) * SDControl^2)/(NThree6 + NControl - 2), 
             J = 1 - 3/(4*dfev - 1), 
             d = (MeanThree6 - MeanControl)/sqrt(pooled_var), 
             g = d * J, 
             vd = (NThree6 + NControl) / (NThree6 * NControl) + d^2/(2 * (NThree6 + NControl)), 
             vg = vd * J^2, 
             log.oddsratio = d*(pi/sqrt(3)),
             var.log.oddsratio = vd*(pi^2)/3,
             a = ((NThree6+NControl)^2)/(NThree6*NControl),
             r = d/sqrt(d^2+a),
             var.r = ((a^2)*vd)/(((d^2)+a)^3)) 

ic <- mutate(.data = ic,
             pooled_var = ((NContact - 1) * SDContact^2 + (NControl - 1) * SDControl^2)/(NContact + NControl - 2), 
             J = 1 - 3/(4*dfev - 1), 
             d = (MeanContact - MeanControl)/sqrt(pooled_var), 
             g = d * J, 
             vd = (NContact + NControl) / (NContact * NControl) + d^2/(2 * (NContact + NControl)), 
             vg = vd * J^2, 
             log.oddsratio = d*(pi/sqrt(3)),
             var.log.oddsratio = vd*(pi^2)/3,
             a = ((NContact+NControl)^2)/(NContact*NControl),
             r = d/sqrt(d^2+a),
             var.r = ((a^2)*vd)/(((d^2)+a)^3)) 


mag <- mutate(.data = mag,
             pooled_var = ((NFemale - 1) * SDFemale^2 + (NControl - 1) * SDControl^2)/(NFemale + NControl - 2), 
             J = 1 - 3/(4*dfev - 1), 
             d = (MeanFemale - MeanControl)/sqrt(pooled_var), 
             g = d * J, 
             vd = (NFemale + NControl) / (NFemale * NControl) + d^2/(2 * (NFemale + NControl)), 
             vg = vd * J^2, 
             log.oddsratio = d*(pi/sqrt(3)),
             var.log.oddsratio = vd*(pi^2)/3,
             a = ((NFemale+NControl)^2)/(NFemale*NControl),
             r = d/sqrt(d^2+a),
             var.r = ((a^2)*vd)/(((d^2)+a)^3)) 

mp <- mutate(.data = mp,
             pooled_var = ((NMoney - 1) * SDMoney^2 + (NControl - 1) * SDControl^2)/(NMoney + NControl - 2), 
             J = 1 - 3/(4*dfev - 1), 
             d = (MeanMoney - MeanControl)/sqrt(pooled_var), 
             g = d * J, 
             vd = (NMoney + NControl) / (NMoney * NControl) + d^2/(2 * (NMoney + NControl)), 
             vg = vd * J^2, 
             log.oddsratio = d*(pi/sqrt(3)),
             var.log.oddsratio = vd*(pi^2)/3,
             a = ((NMoney+NControl)^2)/(NMoney*NControl),
             r = d/sqrt(d^2+a),
             var.r = ((a^2)*vd)/(((d^2)+a)^3)) 

qa <- mutate(.data = qa,
             pooled_var = ((NLiked - 1) * SDLiked^2 + (NControl - 1) * SDControl^2)/(NLiked + NControl - 2), 
             J = 1 - 3/(4*dfev - 1), 
             d = (MeanLiked - MeanControl)/sqrt(pooled_var), 
             g = d * J, 
             vd = (NLiked + NControl) / (NLiked * NControl) + d^2/(2 * (NLiked + NControl)), 
             vg = vd * J^2, 
             log.oddsratio = d*(pi/sqrt(3)),
             var.log.oddsratio = vd*(pi^2)/3,
             a = ((NLiked+NControl)^2)/(NLiked*NControl),
             r = d/sqrt(d^2+a),
             var.r = ((a^2)*vd)/(((d^2)+a)^3)) 

sc <- mutate(.data = sc,
             pooled_var = ((NPaid - 1) * SDPaid^2 + (NFree - 1) * SDFree^2)/(NPaid + NFree - 2), 
             J = 1 - 3/(4*dfev - 1), 
             d = (MeanPaid - MeanFree)/sqrt(pooled_var), 
             g = d * J, 
             vd = (NPaid + NFree) / (NPaid * NFree) + d^2/(2 * (NPaid + NFree)), 
             vg = vd * J^2, 
             log.oddsratio = d*(pi/sqrt(3)),
             var.log.oddsratio = vd*(pi^2)/3,
             a = ((NPaid+NFree)^2)/(NPaid*NFree),
             r = d/sqrt(d^2+a),
             var.r = ((a^2)*vd)/(((d^2)+a)^3)) 

################################
###Coverting from Correlation###
################################

iat <- mutate(.data = iat,
              var.r = ((1-r^2)^2)/(N-1), #Variance of correlation
              z = .5*log((1+r)/(1-r)), #Fisher's Z
              var.z = 1/(N-3), #Variance of Z
              d = 2*r/sqrt(1-r^2), #Cohen's D
              vd = 4*var.r/((1-r^2)^3), #Variance of D
              g = d*(1 - 3/(4*(N - 2) - 1)), #Hedges' G
              vg = vd*((1 - 3/(4*(N - 2) - 1))^2), #Variance of G
              log.oddsratio = d*(pi/sqrt(3)), #Log odds ratio
              var.log.oddsratio = vd*(pi^2)/3) #Variance of log odds ratio

#########################################
###Selecting standardized effect sizes###
#########################################

af <- select(af, Experiment, Site, es.measurement, d, vd, g, vg, r, var.r)
a1 <- select(a1, Experiment, Site, es.measurement, d, vd, g, vg, r, var.r)
a2 <- select(a2, Experiment, Site, es.measurement, d, vd, g, vg, r, var.r)
a3 <- select(a3, Experiment, Site, es.measurement, d, vd, g, vg, r, var.r)
a4 <- select(a4, Experiment, Site, es.measurement, d, vd, g, vg, r, var.r)
fp <- select(fp, Experiment, Site, es.measurement, d, vd, g, vg, r, var.r)
gl <- select(gl, Experiment, Site, es.measurement, d, vd, g, vg, r, var.r)
gf <- select(gf, Experiment, Site, es.measurement, d, vd, g, vg, r, var.r)
iat <- select(iat, Experiment, Site, es.measurement, d, vd, g, vg, r, var.r)
ic <- select(ic, Experiment, Site, es.measurement, d, vd, g, vg, r, var.r)
mag <- select(mag, Experiment, Site, es.measurement, d, vd, g, vg, r, var.r)
mp <- select(mp, Experiment, Site, es.measurement, d, vd, g, vg, r, var.r)
qa <- select(qa, Experiment, Site, es.measurement, d, vd, g, vg, r, var.r)
r <- select(r, Experiment, Site, es.measurement, d, vd, g, vg, r, var.r)
s  <- select(s, Experiment, Site, es.measurement, d, vd, g, vg, r, var.r)
sc <- select(sc, Experiment, Site, es.measurement, d, vd, g, vg, r, var.r)

df <- bind_rows(af, a1, a2, a3, a4, fp, gl, gf, iat, ic, mag, mp, qa, r, s, sc)

names(df)[1:3] = c('experiment', 'site', 'es')
names(df)[ncol(df)] = "vr"

#Writing CSV file
write.csv(df,  "../../manylabs_replicates.csv", row.names=F)
