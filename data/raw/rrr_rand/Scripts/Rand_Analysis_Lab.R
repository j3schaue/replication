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
  data <- data[which(Age>=18 & Age<=34),] # exclude all outside 18-34 and no age
  return(data)
}

# Function for converting a column of Qualtrix output to useable form #
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
                  Moderator=data.frame(D=b1,SE=se.b1))
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
    } else out <- ttest(Contribution,Condition)
    
  }

  # Moderator Analysis #
  if(type=='moderator'){
    
    if(nrow(data)<=4){
      cat('\n*Not enough data! Analyses not feasible.')
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
    
    }
  }
}

file.csv <- list.files(pattern="*.csv")[1] # name of csv file in source folder
# if there are multiple csv files in folder, use first one according to alphabetical order
author <- strsplit(file.csv,split='_',fixed=TRUE)[[1]][1] # extract author's name
data <- read.csv(file.csv,header=FALSE) # read csv file
if(as.character(data[1,1])=='startDate') data <- data[,c(1:10,14:22,11,23:27,12,28:30,13,31:33,13,34:87)] # make new data format compatible
Exclusions <- c('None','(1) Experienced','(2) Non-Compliant','(3) Non-Comprehending','(1) or (2) or (3)')

out <- paste(author,'_Main.txt',sep='') # start writing named output file
sink(out)
cat('MAIN ANALYSIS\n\n')
cat('Outcome Variable\n')
cat('----------------\n')
cat('Contribution: 0-100%\n\n')
cat('Main Predictor\n')
cat('--------------\n')
cat('Condition: TP = Time Pressure, FD = Forced Delay\n\n')
cat('Exclusion Criteria\n')
cat('------------------\n')
cat('(1) Experienced: Participation in Similar Experiments\n')
cat('(2) Non-Compliant: Time >= 10s for TP, Time < 10s for FD\n')
cat('(3) Non-Comprehending: Un1 < 9 or Un2 > 1\n')
cat('\n______________________________________________________\n')
for(i in 1:5){
  cat('\nExclusion: ',Exclusions[i],'\n',rep('=',nchar(Exclusions[i])+11),'\n\n',sep='')
  Rand.Analysis(data,type='main',exclude=i-1)
}
sink()

for(i in c(1,5)){
  if(i==1) Excl <- 'None'
  if(i==5) Excl <- 'All'
  out <- paste(author,'_Moderator_Exclude_',Excl,'.txt',sep='') # start writing named output file
  sink(out)
  cat('MODERATOR ANALYSIS\n\n')
  cat('Outcome Variable\n')
  cat('----------------\n')
  cat('Contribution: 0 to 100%\n\n')
  cat('Main Predictor\n')
  cat('--------------\n')
  cat('Condition: Time Pressure (TP), Forced Delay (FD)\n\n')
  cat('Moderators\n')
  cat('----------\n')
  cat('*Trust: 1 to 10\n')
  cat('Gender: Male, Female, Other\n')
  cat('*Age: (Current Year)-(Birth Year)\n')
  cat('Subject: Number of Subject Pool Experiments\n')
  cat('Paid: Number of Paid Experiments\n')
  cat('MTurk: Number of MTurk Experiments\n')
  cat('Deception: Yes, No to Participation in Deception Experiments\n')
  cat('Knowledge of Other Participants (KOP): Yes, No\n')
  cat('*Horizontal-Individual (HI): 1 to 9\n')
  cat('*Vertical-Individual (VI): 1 to 9\n')
  cat('*Horizontal-Collective (HC): 1 to 9\n')
  cat('*Vertical-Collective (VC): 1 to 9 \n\n')
  cat('*Continuous moderators are mean-centered for analyses\n')
  cat(' (except for Randomization Check)\n\n')
  if(i==1) cat('No Exclusions\n\n')
  if(i==5) {
    cat('All Exclusions\n')
    cat('--------------\n')
    cat('(1) Experienced: Participation in Similar Experiments\n')
    cat('(2) Non-Compliant: Time >= 10s for TP, Time < 10s for FD\n')
    cat('(3) Non-Comprehending: Un1 < 9 or Un2 > 1\n\n')
  }
  cat('________________________________________________________________\n')
  Rand.Analysis(data,type='moderator',exclude=i-1)
  sink()
}
