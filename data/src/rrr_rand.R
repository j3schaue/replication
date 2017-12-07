setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../raw/rrr_rand/Data/")

#Loading required library
library(dplyr)

###################################
## Fuctions from Rand_Analysis.R ##
###################################

# Function for converting a column of Qualtrics output to useable form
convert <- function(x) as.numeric(as.character(x)) 

# Function for removing subjects with missing contribution
rm.missing <- function(data){ 
  TP.missing <- which((is.na(convert(data[,18]))==is.na(convert(data[,19])))==FALSE) # subjects with missing contribution for TP
  FD.missing <- which((is.na(convert(data[,24]))==is.na(convert(data[,25])))==FALSE) # subjects with missing contribution for FD
  missing <- c(TP.missing,FD.missing) # all subjects with missing contribution
  if(length(missing)>0){data <- data[-missing,]} # remove all subjects with missing contribution
  return(data)
}

# Function for removing subjects outside of age 18-34
rm.age <- function(data){
  Year.Born <- convert(data[,85]) + 1919 # birth year
  Age <- 2015 - Year.Born # estimated age
  data <- data[which(Age>=18 & Age<=34),] # exclude all outside 18-34 and no age
  return(data)
}

# Function for t-test 
ttest <- function(C,G){ #C = Contribution, G = Group
  model <- lm(C~G) # C = b0 + b1*G
  df <- model$df
  results <- summary(model)$coef
  b0 <- results[1,1]
  b1 <- results[2,1]
  se.b1 <- results[2,2]
  t.b1 <- results[2,3]
  p.b1 <- results[2,4]
  D <- 2*b1
  SE <- 2*se.b1
  cn <- colnames(results)
  cn <- c(cn[1:3],'df',cn[4])
  test <- round(data.frame(D,SE,t.b1,df,p.b1),4)
  dimnames(test) <- list('',cn)
  TP <- b0+b1
  FD <- b0-b1
  means <- round(data.frame(TP,FD,row.names=''),4)
  CI <- round(data.frame(Lower=D-qt(.975,df)*SE,Upper=D+qt(.975,df)*SE,row.names=''),4)
  out <-  data.frame(TP,FD, D, SE, df, t.b1)
}

# Program
Rand.Analysis.mod <- function(data,type='main',exclude=0) {
  
  #=================#
  # Data Formatting #
  #=================#
  
  # Data #
  colnames(data) <- sapply(data[1,],as.character) # label columns
  data <- data[-(1:2),] # delete first two irrelevant rows
  data <- rm.missing(data) # remove any subjects with missing contribution
  data <- rm.age(data) # remove any subjects outside of age 18-34
  
  # Exclusion Criteria #
  Exp <- convert(data[,67]) # 1 = naive
  TP.Time <- convert(data[,22]) # <10 = compliant
  TP.Time <- ifelse(TP.Time==0,NA,TP.Time) # 0 = NA
  FD.Time <- convert(data[,28]) # >=10 = compliant
  FD.Time <- ifelse(FD.Time==0,NA,FD.Time) # 0 = NA
  Comp1 <- convert(data[,31]) # 9 = comprehend
  Comp2 <- convert(data[,32]) # 1 = comprehend
  if(exclude==1) data <- data[which(Exp==1),] # data excluding experienced
  if(exclude==2) data <- data[which(TP.Time<10 | FD.Time>=10),] # data excluding non-compliant
  if(exclude==3) data <- data[which(Comp1==9 & Comp2==1),] # data excluding non-comprehending
  if(exclude==4) data <- data[which(Exp==1 & (TP.Time<10 | FD.Time>=10) & Comp1==9 & Comp2==1),] # data excluding all
  
  #==========#
  # Analyses #
  #==========#
  
  # Main Variables #
  TP0 <- convert(data[,18]) # $ contribution for time pressure condition
  FD0 <- convert(data[,24]) # $ contribution for forced delay condition
  total <- max(c(TP0,FD0),na.rm=TRUE) # maximum contribution
  TP <- TP0/total*100 # % contribution for time pressure condition
  FD <- FD0/total*100 # % contribution for forced delay condition
  Contribution <- rowSums(cbind(TP,FD),na.rm=TRUE) # TP and FD merged
  Condition <- ifelse(!is.na(TP),'TP','FD') # create Condition vector
  Condition <- factor(Condition,levels=c('TP','FD')) # designate Condition as factor object
  contrasts(Condition) <- contr.sum(2) # TP = 1, FD = -1
  Trust <- convert(data[,73]) # level of trust
  Result <- data.frame(
    nt=sum(Condition=="TP")+1, nf=sum(Condition=="FD"),
    yt=mean(Contribution[which(Condition=="TP")]),
    yf=mean(Contribution[which(Condition=="FD")]),
    vt=var(Contribution[which(Condition=="TP")]),
    vf=var(Contribution[which(Condition=="FD")]))
   
  return(Result)
  }


#############################################################
## Generating ES from raw data with Rand.Analysis function ##
#############################################################

temp = list.files(pattern="*.csv")
for (i in 1:length(temp)) assign(temp[i], read.csv(temp[i], fileEncoding = "latin1"))
original <- read.table("Original_data.txt", header = TRUE)

ac <- Rand.Analysis.mod(Aczel_PAEX_Rand_Data.csv)
ac <- mutate(ac, site = "Aczel")
be <- Rand.Analysis.mod(Begue_PAEX_Rand_Data.csv)
be <- mutate(be, site = "Bègue")
bo <- Rand.Analysis.mod(Bouwmeester_PAEX_Rand_Data.csv)
bo <- mutate(bo, site = "Bouwmeester")
ep <- Rand.Analysis.mod(Espin_PAEX_Rand_Data.csv)
ep <- mutate(ep, site = "Espin")
ev <- Rand.Analysis.mod(Evans_PAEX_Rand_Data.csv)
ev <- mutate(ev, site = "Evans")
fe <- Rand.Analysis.mod(`Ferreira-Santos_PAEX_Rand_Data.csv`)
fe <- mutate(fe, site = "Ferreira-Santos")
fi <- Rand.Analysis.mod(Fiedler_PAEX_Rand_Data.csv)
fi <- mutate(fi, site = "Fiedler")
ha <- Rand.Analysis.mod(Hauser_PAEX_Rand_Data.csv)
ha <- mutate(ha, site = "Hauser")
he <- Rand.Analysis.mod(Hernan_PAEX_Rand_Data.csv)
he <- mutate(he, site = "Hernan")
lo <- Rand.Analysis.mod(Lohse_PAEX_Rand_Data.csv)
lo <- mutate(lo, site = "Lohse")
mi <- Rand.Analysis.mod(Mischkowski_PAEX_Rand_Data.csv)
mi <- mutate(mi, site = "Mischkowski")
ne <- Rand.Analysis.mod(Neal_PAEX_Rand_Data.csv)
ne <- mutate(ne, site = "Neal")
no <- Rand.Analysis.mod(Novakva_PAEX_Rand_Data.csv)
no <- mutate(no, site = "Novakva")
pa <- Rand.Analysis.mod(Paga_PAEX_Rand_Data.csv)
pa <- mutate(pa, site = "Paga")
pi <- Rand.Analysis.mod(Piovesan_PAEX_Rand_Data.csv)
pi <- mutate(pi, site = "Piovesan")
sa <- Rand.Analysis.mod(Salomon_PAEX_Rand_Data.csv)
sa <- mutate(sa, site = "Salomon")
sr <- Rand.Analysis.mod(Srinivasan_PAEX_Rand_Data.csv)
sr <- mutate(sr, site = "Srinivasan")
ri <- Rand.Analysis.mod(Tinghog_PAEX_Rand_Data.csv)
ri <- mutate(ri, site = "Tinghog")
tr <- Rand.Analysis.mod(Trueblood_PAEX_Rand_Data.csv)
tr <- mutate(tr, site = "Trueblood")
wi <- Rand.Analysis.mod(Wills_PAEX_Rand_Data.csv)
wi <- mutate(wi, site = "Wills")
wo <- Rand.Analysis.mod(Wollbrant_PAEX_Rand_Data.csv)
wo <- mutate(wo, site = "Wollbrant")

# original study
oSE=5.316327 # from RRR paper
odiff = 8.58
ont=55; onf=98 # sample sizes from Rand (2012) supplement p. 15
od = odiff/sqrt(oSE^2 * ont * onf/(ont + onf))
orig = data.frame(nt=ont, nf=98, yt=NA, yf=NA, vt=NA, vf=NA, site='original',
                  d=od,
                  vd=(ont + onf)/(ont*onf) + od^2/(2*(onf + ont)), 
                  es='smd', 
                  n = ont + 98)

# from collin. excluded
# # original study from Rand, Greene, and Novak (2012)
# or = data.frame(D=8.58, SE=5.316327, TP=)
# # or <- c(D=49.43, 40.85, -8.6, 5.274374, 153, -1.63, "original") #from Rand paper

# Full data frame
rrr_rand <- rbind.data.frame(ac, be, bo, ep, ev, fe, fi, ha, he, lo, mi, ne, no, pa, pi, sa, sr, ri, tr, wi, wo) %>%
  mutate(d = (yt - yf)/sqrt((nt*vt + nf*vf)/(nf + nt - 2)), 
         vd = (nt + nf)/(nt * nf) + d^2/(2*(nt + nf)), 
         es = 'smd', 
         n = nt + nf) %>% 
  rbind(., orig)

rrr_rand$experiment = "time/delay"
rrr_rand$replicated = 0

tail(rrr_rand)

write.csv(rrr_rand, "../../../rrr_rand.csv", row.names = F)
