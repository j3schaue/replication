#################
# RAND ANALYSES #
#################

# Function for removing subjects with missing contribution #
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
  data <- data[which(Age>=18 & Age<=34),] # exclude all outside 18-34
  return(data)
}

# Function for converting a column of Qualtrics output to useable form #
convert <- function(x) as.numeric(as.character(x)) 

# Function for t-test #
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
  cat('Mean Contribution per Condition\n')
  cat('-------------------------------\n')
  print(means)
  cat('\nt-test: Difference in Means\n')
  cat('---------------------------\n')
  print(test)
  cat('\n95% CI: Difference in Means\n')
  cat('---------------------------\n')
  print(CI)
  cat('\n______________________________________________________\n')
  
  out <-  data.frame(D,SE,TP,FD)
  
}

# Function for ANCOVA #
ancova <- function(C,G,X){ #C = Contribution, G = Group, X = Moderator
  mod <- deparse(substitute(X))
  
  if(is.numeric(X)){
    Xm.TP <- mean(X[G=='TP'],na.rm=TRUE)
    Xm.FD <- mean(X[G=='FD'],na.rm=TRUE)
    mod.m <- round(data.frame(Xm.TP,Xm.FD),4)
    dimnames(mod.m) <- list('',c('TP','FD'))
  } else{
    if(mod=='Gender'){
      XpM.TP <- mean(X[G=='TP']==1,na.rm=TRUE)
      XpF.TP <- mean(X[G=='TP']==2,na.rm=TRUE)
      XpM.FD <- mean(X[G=='FD']==1,na.rm=TRUE)
      XpF.FD <- mean(X[G=='FD']==2,na.rm=TRUE)
      mod.p <- round(data.frame(matrix(c(XpM.TP,XpF.TP,XpM.FD,XpF.FD),2,2)),4)
      dimnames(mod.p) <- list(c('Male','Female'),c('TP','FD'))
    } else{
      Xp1.TP <- mean(X[G=='TP']==1,na.rm=TRUE)
      Xp0.TP <- mean(X[G=='TP']==0,na.rm=TRUE)
      Xp1.FD <- mean(X[G=='FD']==1,na.rm=TRUE)
      Xp0.FD <- mean(X[G=='FD']==0,na.rm=TRUE)
      mod.p <- round(data.frame(matrix(c(Xp1.TP,Xp0.TP,Xp1.FD,Xp0.FD),2,2)),4)
      dimnames(mod.p) <- list(c('Yes','No'),c('TP','FD'))
    }
  }
  
  fail <- var(as.numeric(X[G=='TP']),na.rm=TRUE) %in% c(0,NA) | var(as.numeric(X[G=='FD']),na.rm=TRUE) %in% c(0,NA)
  
  if(!fail){
    if(is.numeric(X)){
      M <- X-mean(X,na.rm=TRUE) # center continuous moderator
    } else M <- X
    modelA <- lm(C~M*G) # C = b0' + b1'*M + b2'*G + b3*M*G
    dfA <- modelA$df
    resultA <- summary(modelA)$coef
    cn <- colnames(resultA)
    cn <- c(cn[1:3],'df',cn[4])
    b3 <- resultA[4,1]
    se.b3 <- resultA[4,2]
    t.b3 <- resultA[4,3]
    p.b3 <- resultA[4,4]
    b1 <- resultA[2,1]
    b2 <- resultA[3,1]
    if(is.numeric(M)){
      D3 <- 2*b3
      SE3 <- 2*se.b3
      S1 <- b1+b3
      S2 <- b1-b3
    } else{
      D3 <- 4*b3
      SE3 <- 4*se.b3
      S1 <- 2*(b2+b3)
      S2 <- 2*(b2-b3)
    }
    slopes <- round(data.frame(S1,S2,row.names=''),4)
    test3 <- round(data.frame(D3,SE3,t.b3,dfA,p.b3),4)
    dimnames(test3) <- list('',cn)
    CI3 <- round(data.frame(Lower=D3-qt(.975,dfA)*SE3,Upper=D3+qt(.975,dfA)*SE3,row.names=''),4)
    
    modelB <- lm(C~M+G) # C = b0 + b1*M + b2*G
    dfB <- modelB$df
    resultB <- summary(modelB)$coef 
    b0 <- resultB[1,1]
    b1 <- resultB[2,1]
    se.b1 <- resultB[2,2]
    t.b1 <- resultB[2,3]
    p.b1 <- resultB[2,4]
    b2 <- resultB[3,1]
    se.b2 <- resultB[3,2]
    t.b2 <- resultB[3,3]
    p.b2 <- resultB[3,4]
    D2 <- 2*b2
    SE2 <- 2*se.b2
    TP <- b0+b2
    FD <- b0-b2
    means2 <- round(data.frame(TP,FD,row.names=''),4)
    test2 <- round(data.frame(D2,SE2,t.b2,dfB,p.b2),4)
    CI2 <- round(data.frame(Lower=D2-qt(.975,dfB)*SE2,Upper=D2+qt(.975,dfB)*SE2,row.names=''),4)
    if(is.numeric(M)){
      colnames(slopes) <- c('TP','FD')
      test1 <- round(data.frame(b1,se.b1,t.b1,dfB,p.b1),4)
      CI1 <- round(data.frame(Lower=b1-qt(.975,dfB)*se.b1,Upper=b1+qt(.975,dfB)*se.b1,row.names=''),4)
    } else{
      D1 <- 2*b1
      SE1 <- 2*se.b1
      if(mod=='Gender'){
        colnames(slopes) <- c('Male','Female')
        Male <- b0+b1
        Female <- b0-b1
        means1 <- round(data.frame(Male,Female,row.names=''),4)
      } else{
        colnames(slopes) <- c('Yes','No')
        Yes <- b0+b1
        No <- b0-b1
        means1 <- round(data.frame(Yes,No,row.names=''),4)
      }
      test1 <- round(data.frame(D1,SE1,t.b1,dfB,p.b1),4)
      CI1 <- round(data.frame(Lower=D1-qt(.975,dfB)*SE1,Upper=D1+qt(.975,dfB)*SE1,row.names=''),4)
    }
    dimnames(test2) <- dimnames(test1) <- list('',cn)
  }
    
  if(is.numeric(X)){
    cat('\n[',toupper(mod),']\n\n',sep='')
    cat('Randomization Check\n')
    cat('===================\n\n')
    cat('Mean of',mod,'per Condition\n')
    cat(rep('-',nchar(mod)+22),'\n',sep='')
    print(mod.m)
    if(!fail){
      cat('\nModeration Effect\n')
      cat('=================\n\n')
      cat('Model: Contribution ~ ',mod,' + Condition + ',mod,':Condition\n\n',sep='')
      cat('Slope of',mod,'per Condition\n')
      cat(rep('-',nchar(mod)+23),'\n',sep='')
      print(slopes)
      cat('\nt-test: Difference in Slopes\n')
      cat('----------------------------\n')
      print(test3)
      cat('\n95% CI: Difference in Slopes\n')
      cat('----------------------------\n')
      print(CI3)
      cat('\n')
      if(p.b3<.05) {cat('*Warning message:\n Significant moderation effect of',mod,'on effect size.\n Interpret ANCOVA with caution!\n\n')}
      cat('ANCOVA\n')
      cat('======\n\n')
      cat('Model: Contribution ~',mod,'+ Condition\n\n')
      cat('Adjusted Mean per Condition\n')
      cat('---------------------------\n')
      print(means2)
      cat('\nt-test: Difference in Adjusted Means\n')
      cat('------------------------------------\n')
      print(test2)
      cat('\n95% CI: Difference in Adjusted Means\n')
      cat('------------------------------------\n')
      print(CI2)
      cat('\n.........................................\n')
      cat('\nt-test: Slope of',mod,'\n')
      cat(rep('-',nchar(mod)+17),'\n',sep='')
      print(test1)
      cat('\n95% CI: Slope of',mod,'\n')
      cat(rep('-',nchar(mod)+17),'\n',sep='')
      print(CI1)
    } else cat('\n*Randomization fail! Moderator analysis not feasible.\n')
    cat('\n________________________________________________________________\n')
  } else{
    cat('\n[',toupper(mod),']\n\n',sep='')
    cat('Randomization Check\n')
    cat('===================\n\n')
    cat('Proportions of',mod,'per Condition\n')
    cat(rep('-',nchar(mod)+29),'\n',sep='')
    print(mod.p)
    if(!fail){
      if(mod=='Gender') cat('\n*Other gender omitted from analysis\n')
      cat('\nModeration Effect\n')
      cat('=================\n\n')
      cat('Model: Contribution ~ ',mod,' + Condition + ',mod,':Condition\n\n',sep='')
      cat('Effect Size per Level of',mod,'\n')
      cat(rep('-',nchar(mod)+25),'\n',sep='')
      print(slopes)
      cat('\nt-test: Difference in Effect Sizes\n')
      cat('----------------------------------\n')
      print(test3)
      cat('\n95% CI: Difference in Effect Sizes\n')
      cat('----------------------------------\n')
      print(CI3)
      cat('\n')
      if(p.b3<.05) {cat('*Warning message:\n Significant moderation effect of',mod,'on effect size.\n Interpret factor effects with caution!\n\n')}
      cat('Factor Effects (Additive ANOVA)\n')
      cat('===============================\n\n')
      cat('Model: Contribution ~',mod,'+ Condition\n\n')
      cat('Marginal Mean per Condition\n')
      cat('---------------------------\n')
      print(means2)
      cat('\nt-test: Difference in Marginal Means\n')
      cat('------------------------------------\n')
      print(test2)
      cat('\n95% CI: Difference in Marginal Means\n')
      cat('------------------------------------\n')
      print(CI2)
      cat('\n...................................................\n')
      cat('\nMarginal Mean per',mod,'\n')
      cat(rep('-',nchar(mod)+18),'\n',sep='')
      print(means1)
      cat('\nt-test: Difference in Marginal Means\n')
      cat('------------------------------------\n')
      print(test2)
      cat('\n95% CI: Difference in Marginal Means\n')
      cat('------------------------------------\n')
      print(CI1)
    } else cat('\n*Randomization fail! Moderator analysis not feasible.\n')
    cat('\n________________________________________________________________\n')
  } 
  
  if(!fail){
    if(is.numeric(M)) {
      out <- list(Interaction=data.frame(D=D3,SE=SE3,TP=S1,FD=S2),
                  Condition=data.frame(D=D2,SE=SE2,TP,FD),
                  Moderator=data.frame(D=b1,SE=se.b1,M1=NA,M2=NA))
    }
    if(is.factor(M) & mod=='Gender') {
      out <- list(Interaction=data.frame(D=D3,SE=SE3,TP=S1,FD=S2),
                  Condition=data.frame(D=D2,SE=SE2,TP,FD),
                  Moderator=data.frame(D=D1,SE=SE1,Male,Female))
    }
    if(is.factor(M) & mod!='Gender') {
      out <- list(Interaction=data.frame(D=D3,SE=SE3,TP=S1,FD=S2),
                  Condition=data.frame(D=D2,SE=SE2,TP,FD),
                  Moderator=data.frame(D=D1,SE=SE1,Yes,No))
    }
  } else out <- list(Interaction=NA,Condition=NA,Moderator=NA)
  
  out <- out
  
}
  
# Program #
Rand.Analysis <- function(data,type='main',exclude=0) {
  
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
                        
  # Main Analysis #
  if(type=='main'){
    
    if(any(c(sum(Condition=='TP'),sum(Condition=='FD'))<2)){
      cat('*Not enough data! Analyses not feasible.\n')
      cat('\n______________________________________________________\n')
      meta.out <- data.frame(D=NA,SE=NA,TP=NA,FD=NA,Trust=NA)
    } else {
      out <- ttest(Contribution,Condition)
      Trust.m <- mean(Trust,na.rm=TRUE)
      meta.out <- data.frame(out,Trust=Trust.m)
    }

  }

  # Moderator Analysis #
  if(type=='moderator'){
    
    if(nrow(data)<=4){
      cat('\n*Not enough data! Analyses not feasible.')
      out <- list(Interaction=NA,Condition=NA,Moderator=NA)
      meta.out <- list(Trust=out,Gender=out,Age=out,Subject=out,
                       Paid=out,MTurk=out,Deception=out,KOP=out,
                       HI=out,VI=out,HC=out,VC=out)
    } else{
      
      Trust.out <- ancova(Contribution,Condition,Trust)
      
      Gender <- convert(data[,86]) # 1 = male, 2 = female, 3 = other
      Gender[Gender==3] <- NA # remove other
      Gender <- factor(Gender) # designate Gender as factor
      if(!(var(as.numeric(Gender),na.rm=TRUE) %in% c(0,NA))) contrasts(Gender) <- contr.sum(2) # 1 = male, -1 = female
      Gender.out <- ancova(Contribution,Condition,Gender)
      
      #Year.Now <- as.numeric(format(Sys.Date(),format='%Y')) # current year
      Year.Born <- convert(data[,85]) + 1919 # birth year
      Age <- 2015 - Year.Born # estimated age
      Age.out <- ancova(Contribution,Condition,Age)
      
      Subject <- convert(data[,68]) # number of subject pool experiments
      Subject[Subject>0] <- 1 # replace all non-zero values with 1 (0 = no, 1 = yes)
      Subject <- factor(Subject,levels=c(1,0)) # designate Subject as factor
      if(!(var(as.numeric(Subject),na.rm=TRUE) %in% c(0,NA))) contrasts(Subject) <- contr.sum(2) # 1 = male, -1 = female
      Subject.out <- ancova(Contribution,Condition,Subject)
      
      Paid <- convert(data[,69]) # number of paid experiments
      Paid[Paid>0] <- 1 # replace all non-zero values with 1 (0 = no, 1 = yes)
      Paid <- factor(Paid,levels=c(1,0)) # designate Paid as factor
      if(!(var(as.numeric(Paid),na.rm=TRUE) %in% c(0,NA))) contrasts(Paid) <- contr.sum(2) # 1 = yes, -1 = no
      Paid.out <- ancova(Contribution,Condition,Paid)
      
      MTurk <- convert(data[,70]) # number of MTurk experiments
      MTurk[MTurk>0] <- 1 # replace all non-zero values with 1 (0 = no, 1 = yes)
      MTurk <- factor(MTurk,levels=c(1,0)) # designate MTurk as factor
      if(!(var(as.numeric(MTurk),na.rm=TRUE) %in% c(0,NA))) contrasts(MTurk) <- contr.sum(2) # 1 = yes, -1 = no
      MTurk.out <- ancova(Contribution,Condition,MTurk)
      
      Deception <- convert(data[,71]) # number of deception experiments
      Deception[Deception>0] <- 1 # replace all non-zero values with 1 (0 = no, 1 = yes)
      Deception <- factor(Deception,levels=c(1,0)) # designate KOP as factor
      if(!(var(as.numeric(Deception),na.rm=TRUE) %in% c(0,NA))) contrasts(Deception) <- contr.sum(2) # 1 = yes, -1 = no
      Deception.out <- ancova(Contribution,Condition,Deception)
      
      KOP <- convert(data[,88]) # number of participants known
      KOP[KOP>0] <- 1 # replace all non-zero values with 1 (0 = no, 1 = yes)
      KOP <- factor(KOP,levels=c(1,0)) # designate KOP as factor
      if(!(var(as.numeric(KOP),na.rm=TRUE) %in% c(0,NA))) contrasts(KOP) <- contr.sum(2) # 1 = yes, -1 = no
      KOP.out <- ancova(Contribution,Condition,KOP)
      
      HI.i <- sapply(data[,35:42],convert) # item responses on horizontal-individual scale
      HI <- rowMeans(HI.i,na.rm=TRUE) # mean score on horizontal-individual scale
      HI.out <- ancova(Contribution,Condition,HI)
      
      VI.i <- sapply(data[,43:50],convert) # item responses on vertical-individual scale
      VI.i[,8] <- 10-VI.i[,8] # reverse scoring of item 8
      VI <- rowMeans(VI.i,na.rm=TRUE) # mean score on vertical-individual scale
      VI.out <- ancova(Contribution,Condition,VI)
      
      HC.i <- sapply(data[,51:58],convert) # item responses on horizontal-collective scale
      HC <- rowMeans(HC.i,na.rm=TRUE) # mean score on horizontal-collective scale
      HC.out <- ancova(Contribution,Condition,HC)
      
      VC.i <- sapply(data[,59:66],convert) # item responses on vertical-collective scale
      VC <- rowMeans(VC.i,na.rm=TRUE) # mean score on vertical-collective scale
      VC.out <- ancova(Contribution,Condition,VC)
    
      meta.out <- list(Trust=Trust.out,Gender=Gender.out,Age=Age.out,Subject=Subject.out,
                       Paid=Paid.out,MTurk=MTurk.out,Deception=Deception.out,KOP=KOP.out,
                       HI=HI.out,VI=VI.out,HC=HC.out,VC=VC.out)
    }
    
  }
  
  meta.out <- meta.out
  
}

# Descriptive Stats #
Rand.DS <- function(data) {
    
  # Data #
  colnames(data) <- sapply(data[1,],as.character) # label columns
  data <- data[-(1:2),] # delete first two irrelevant rows
  data <- rm.missing(data) # remove any subjects with missing contribution
  data <- rm.age(data) # remove any subjects outside of age 18-34
  Exp <- convert(data[,67]) # 1 = naive
  TP.Time <- convert(data[,22]) # <10 = compliant
  TP.Time <- ifelse(TP.Time==0,NA,TP.Time) # 0 = NA
  FD.Time <- convert(data[,28]) # >=10 = compliant
  FD.Time <- ifelse(FD.Time==0,NA,FD.Time) # 0 = NA
  Comp1 <- convert(data[,31]) # 9 = comprehend
  Comp2 <- convert(data[,32]) # 1 = comprehend
  Condition <- ifelse(!is.na(TP.Time),'TP','FD') # create Condition vector

  # Moderators #
  N_Total_TP <- sum(!is.na(TP.Time)) # total sample size of time pressure condition
  N_Total_FD <- sum(!is.na(FD.Time)) # total sample size of forced delay condition
  
  N_Women_TP <- sum(convert(data[,86])[Condition=='TP']==2,na.rm=TRUE) # number of women in time pressure condition
  N_Women_FD <- sum(convert(data[,86])[Condition=='FD']==2,na.rm=TRUE) # number of women in forced delay condition
  
  #Year.Now <- as.numeric(format(Sys.Date(),format='%Y')) # current year
  Year.Born <- convert(data[,85]) + 1919 # birth year
  Age <- 2015 - Year.Born # estimated age
  M_Age_TP <- mean(Age[Condition=='TP'],na.rm=TRUE) # mean age of time pressure condition
  SD_Age_TP <- sd(Age[Condition=='TP'],na.rm=TRUE) # standard deviation of age of time pressure condition
  M_Age_FD <- mean(Age[Condition=='FD'],na.rm=TRUE) # mean age of time pressure condition
  SD_Age_FD <- sd(Age[Condition=='FD'],na.rm=TRUE) # standard deviation of age of time pressure condition    
  
  N_Personal_TP <- sum(Comp2[Condition=='TP']==1,na.rm=TRUE) # number comprehending personal benefit for time pressure condition
  N_Personal_FD <- sum(Comp2[Condition=='FD']==1,na.rm=TRUE) # number comprehending personal benefit for forced delay condition
  N_Group_TP <- sum(Comp1[Condition=='TP']==9,na.rm=TRUE) # number comprehending group benefit for time pressure condition
  N_Group_FD <- sum(Comp1[Condition=='FD']==9,na.rm=TRUE) # number comprehending group benefit for forced delay condition

  M_Trust1_TP <- mean(convert(data[,72])[Condition=='TP'],na.rm=TRUE) # mean Trust1 for time pressure condition
  SD_Trust1_TP <- sd(convert(data[,72])[Condition=='TP'],na.rm=TRUE) # standard deviation of Trust1 for forced delay condition
  M_Trust1_FD <- mean(convert(data[,72])[Condition=='FD'],na.rm=TRUE) # mean Trust1 for time pressure condition
  SD_Trust1_FD <- sd(convert(data[,72])[Condition=='FD'],na.rm=TRUE) # standard deviation of Trust1 for forced delay condition
  
  M_Trust2_TP <- mean(convert(data[,73])[Condition=='TP'],na.rm=TRUE) # mean Trust2 for time pressure condition
  SD_Trust2_TP <- sd(convert(data[,73])[Condition=='TP'],na.rm=TRUE) # standard deviation of Trust2 for forced delay condition
  M_Trust2_FD <- mean(convert(data[,73])[Condition=='FD'],na.rm=TRUE) # mean Trust2 for time pressure condition
  SD_Trust2_FD <- sd(convert(data[,73])[Condition=='FD'],na.rm=TRUE) # standard deviation of Trust2 for forced delay condition
  
  N_Trust3_TP <- sum(convert(data[,74])[Condition=='TP']==2,na.rm=TRUE) # number of low Trust3 in time pressure condition
  N_Trust3_FD <- sum(convert(data[,74])[Condition=='FD']==2,na.rm=TRUE) # number of low Trust3 in forced delay condition
  
  N_Trust4_TP <- sum(convert(data[,75])[Condition=='TP']==1,na.rm=TRUE) # number of low Trust4 in time pressure condition
  N_Trust4_FD <- sum(convert(data[,75])[Condition=='FD']==1,na.rm=TRUE) # number of low Trust4 in forced delay condition
  
  N_Trust5_TP <- sum(convert(data[,76])[Condition=='TP']==2,na.rm=TRUE) # number of low Trust5 in time pressure condition
  N_Trust5_FD <- sum(convert(data[,76])[Condition=='FD']==2,na.rm=TRUE) # number of low Trust5 in forced delay condition
  
  N_Trust6_TP <- sum(convert(data[,77])[Condition=='TP']==1,na.rm=TRUE) # number of low Trust6 in time pressure condition
  N_Trust6_FD <- sum(convert(data[,77])[Condition=='FD']==1,na.rm=TRUE) # number of low Trust6 in forced delay condition
  
  N_Exp1_TP <- sum(convert(data[,67])[Condition=='TP']==1,na.rm=TRUE) # number naive in time pressure condition
  N_Exp1_FD <- sum(convert(data[,67])[Condition=='FD']==1,na.rm=TRUE) # number naive in forced delay condition
  
  M_Exp2_TP <- mean(convert(data[,68])[Condition=='TP'],na.rm=TRUE) # mean subject pool experiments in time pressure condition
  SD_Exp2_TP <- sd(convert(data[,68])[Condition=='TP'],na.rm=TRUE) # standard deviation of subject pool experiments in time pressure condition
  M_Exp2_FD <- mean(convert(data[,68])[Condition=='FD'],na.rm=TRUE) # mean subject pool experiments in forced delay condition
  SD_Exp2_FD <- sd(convert(data[,68])[Condition=='FD'],na.rm=TRUE) # standard deviation of subject pool experiments in forced delay condition
  
  M_Exp3_TP <- mean(convert(data[,69])[Condition=='TP'],na.rm=TRUE) # mean paid experiments in time pressure condition
  SD_Exp3_TP <- sd(convert(data[,69])[Condition=='TP'],na.rm=TRUE) # standard deviation of paid experiments in time pressure condition
  M_Exp3_FD <- mean(convert(data[,69])[Condition=='FD'],na.rm=TRUE) # mean paid experiments in forced delay condition
  SD_Exp3_FD <- sd(convert(data[,69])[Condition=='FD'],na.rm=TRUE) # standard deviation of paid experiments in forced delay condition
  
  M_Exp4_TP <- mean(convert(data[,70])[Condition=='TP'],na.rm=TRUE) # mean MTurk experiments in time pressure condition
  SD_Exp4_TP <- sd(convert(data[,70])[Condition=='TP'],na.rm=TRUE) # standard deviation of MTurk experiments in time pressure condition
  M_Exp4_FD <- mean(convert(data[,70])[Condition=='FD'],na.rm=TRUE) # mean MTurk experiments in forced delay condition
  SD_Exp4_FD <- sd(convert(data[,70])[Condition=='FD'],na.rm=TRUE) # standard deviation of MTurk experiments in forced delay condition
      
  M_Exp5_TP <- mean(convert(data[,71])[Condition=='TP'],na.rm=TRUE) # mean deception experiments in time pressure condition
  SD_Exp5_TP <- sd(convert(data[,71])[Condition=='TP'],na.rm=TRUE) # standard deviation of deception experiments in time pressure condition
  M_Exp5_FD <- mean(convert(data[,71])[Condition=='FD'],na.rm=TRUE) # mean deception experiments in forced delay condition
  SD_Exp5_FD <- sd(convert(data[,71])[Condition=='FD'],na.rm=TRUE) # standard deviation of deception experiments in forced delay condition
  
  M_KOP_TP <- mean(convert(data[,88])[Condition=='TP'],na.rm=TRUE) # mean participants known in time pressure condition
  SD_KOP_TP <- sd(convert(data[,88])[Condition=='TP'],na.rm=TRUE) # standard deviation of participants known in time pressure condition
  M_KOP_FD <- mean(convert(data[,88])[Condition=='FD'],na.rm=TRUE) # mean participants known in forced delay condition
  SD_KOP_FD <- sd(convert(data[,88])[Condition=='FD'],na.rm=TRUE) # standard deviation of participants known in forced delay condition
      
  HI.i <- sapply(data[,35:42],convert) # item responses on horizontal-individual scale
  HI <- rowMeans(HI.i,na.rm=TRUE) # mean score on horizontal-individual scale
  M_HI_TP <- mean(HI[Condition=='TP'],na.rm=TRUE) # mean HI in time pressure condition
  SD_HI_TP <- sd(HI[Condition=='TP'],na.rm=TRUE) # standard deviation of HI in time pressure condition
  M_HI_FD <- mean(HI[Condition=='FD'],na.rm=TRUE) # mean HI in forced delay condition
  SD_HI_FD <- sd(HI[Condition=='FD'],na.rm=TRUE) # standard deviation of HI in forced delay condition

  VI.i <- sapply(data[,43:50],convert) # item responses on vertical-individual scale
  VI.i[,8] <- 10-VI.i[,8] # reverse scoring of item 8
  VI <- rowMeans(VI.i,na.rm=TRUE) # mean score on vertical-individual scale
  M_VI_TP <- mean(VI[Condition=='TP'],na.rm=TRUE) # mean VI in time pressure condition
  SD_VI_TP <- sd(VI[Condition=='TP'],na.rm=TRUE) # standard deviation of VI in time pressure condition
  M_VI_FD <- mean(VI[Condition=='FD'],na.rm=TRUE) # mean VI in forced delay condition
  SD_VI_FD <- sd(VI[Condition=='FD'],na.rm=TRUE) # standard deviation of VI in forced delay condition
      
  HC.i <- sapply(data[,51:58],convert) # item responses on horizontal-collective scale
  HC <- rowMeans(HC.i,na.rm=TRUE) # mean score on horizontal-collective scale
  M_HC_TP <- mean(HC[Condition=='TP'],na.rm=TRUE) # mean HC in time pressure condition
  SD_HC_TP <- sd(HC[Condition=='TP'],na.rm=TRUE) # standard deviation of HC in time pressure condition
  M_HC_FD <- mean(HC[Condition=='FD'],na.rm=TRUE) # mean HC in forced delay condition
  SD_HC_FD <- sd(HC[Condition=='FD'],na.rm=TRUE) # standard deviation of HC in forced delay condition
      
  VC.i <- sapply(data[,59:66],convert) # item responses on vertical-collective scale
  VC <- rowMeans(VC.i,na.rm=TRUE) # mean score on vertical-collective scale
  M_VC_TP <- mean(VC[Condition=='TP'],na.rm=TRUE) # mean VC in time pressure condition
  SD_VC_TP <- sd(VC[Condition=='TP'],na.rm=TRUE) # standard deviation of VC in time pressure condition
  M_VC_FD <- mean(VC[Condition=='FD'],na.rm=TRUE) # mean VC in forced delay condition
  SD_VC_FD <- sd(VC[Condition=='FD'],na.rm=TRUE) # standard deviation of VC in forced delay condition
  
  stats1.out <- c(N_Total_TP,N_Women_TP,M_Age_TP,SD_Age_TP,N_Personal_TP,N_Group_TP,M_Trust1_TP,SD_Trust1_TP,M_Trust2_TP,SD_Trust2_TP,
                  N_Trust3_TP,N_Trust4_TP,N_Trust5_TP,N_Trust6_TP,N_Exp1_TP,M_Exp2_TP,SD_Exp2_TP,M_Exp3_TP,SD_Exp3_TP,M_Exp4_TP,SD_Exp4_TP,
                  M_Exp5_TP,SD_Exp5_TP,M_KOP_TP,SD_KOP_TP,M_HI_TP,SD_HI_TP,M_VI_TP,SD_VI_TP,M_HC_TP,SD_HC_TP,M_VC_TP,SD_VC_TP,
                  N_Total_FD,N_Women_FD,M_Age_FD,SD_Age_FD,N_Personal_FD,N_Group_FD,M_Trust1_FD,SD_Trust1_FD,M_Trust2_FD,SD_Trust2_FD,
                  N_Trust3_FD,N_Trust4_FD,N_Trust5_FD,N_Trust6_FD,N_Exp1_FD,M_Exp2_FD,SD_Exp2_FD,M_Exp3_FD,SD_Exp3_FD,M_Exp4_FD,SD_Exp4_FD,
                  M_Exp5_FD,SD_Exp5_FD,M_KOP_FD,SD_KOP_FD,M_HI_FD,SD_HI_FD,M_VI_FD,SD_VI_FD,M_HC_FD,SD_HC_FD,M_VC_FD,SD_VC_FD)
  
  # Main #
  stats2 <- matrix(0,10,5)
  
  for(exclude in 0:4){
    
    # Exclusion Criteria #  
    if(exclude==0) data.c <- data # data excluding none
    if(exclude==1) data.c <- data[which(Exp==1),] # data excluding experienced
    if(exclude==2) data.c <- data[which(TP.Time<10 | FD.Time>=10),] # data excluding non-compliant
    if(exclude==3) data.c <- data[which(Comp1==9 & Comp2==1),] # data excluding non-comprehending
    if(exclude==4) data.c <- data[which(Exp==1 & (TP.Time<10 | FD.Time>=10) & Comp1==9 & Comp2==1),] # data excluding all
    
    # Contribution #
    TP0 <- convert(data.c[,18]) # $ contribution for time pressure condition
    FD0 <- convert(data.c[,24]) # $ contribution for forced delay condition
    total <- max(c(TP0,FD0),na.rm=TRUE) # maximum contribution
    TP <- TP0/total*100 # % contribution for time pressure condition
    FD <- FD0/total*100 # % contribution for forced delay condition
    TP.T <- convert(data.c[,22]) # decision time for time pressure condition
    TP.T <- ifelse(TP.T==0,NA,TP.T) # 0 = NA
    FD.T <- convert(data.c[,28]) # decision time for forced delay condition
    FD.T <- ifelse(FD.T==0,NA,FD.T) # 0 = NA
    
    # Stats #
    N_TP <- sum(!is.na(TP)) # sample size of time pressure condition
    DT_M_TP <- mean(TP.T,na.rm=TRUE) # mean decision time for time pressure condition
    DT_SD_TP <- sd(TP.T,na.rm=TRUE) # standard deviation of decision time for time pressure condition
    C_M_TP <- mean(TP,na.rm=TRUE) # mean % contribution for time pressure condition 
    C_SD_TP <- sd(TP,na.rm=TRUE) # standard deviation of % contribution for time pressure condition
    N_FD <- sum(!is.na(FD)) # sample size of forced delay condition
    DT_M_FD <- mean(FD.T,na.rm=TRUE) # mean decision time for forced delay condition
    DT_SD_FD <- sd(FD.T,na.rm=TRUE) # standard deviation of decision time for forced delay condition
    C_M_FD <- mean(FD,na.rm=TRUE) # mean % contribution for forced delay condition
    C_SD_FD <- sd(FD,na.rm=TRUE) # standard deviation of % contribution for forced delay condition
    stats2[,exclude+1] <- c(N_TP,DT_M_TP,DT_SD_TP,C_M_TP,C_SD_TP,N_FD,DT_M_FD,DT_SD_FD,C_M_FD,C_SD_FD)
    
  }
  
  stats2.out <- matrix(stats2,1,50)
  
  out <- round(c(stats1.out,stats2.out),1)
      
}
    

