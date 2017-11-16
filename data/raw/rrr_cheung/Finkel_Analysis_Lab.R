###################
# FINKEL ANALYSES #
###################

cat('\nPlease be patient. This may take a few minutes.\n')

convert <- function(x) as.numeric(as.character(x)) # function that converts a column of Qualtrix output to useable form in R

Finkel.Analyze <- function(data,include=0:2,r=10000) {

#=================#
# Data Formatting #
#=================#

# NOTE: Subject's mean score on a scale is considered missing (NA) if any of the items on the scale is unanswered

# Read Dataset #
colnames(data) <- apply(data[1,],2,as.character) # label columns
data <- data[-(1:2),] # delete first two irrelevant rows
Exclude <- convert(data$Exclude) # exclusion criteria
data <- data[which(Exclude %in% include),] # include specified observations

# Gender #
Gender <- with(data,convert(Q5.1)) # convert gender item
Gender[Gender==1] <- 'M' # recode 1 (male) --> M
Gender[Gender==2] <- 'F' # recode 2 (female) --> F
Gender <- factor(Gender,levels=c('M','F'))  # designate Gender as factor object
contrasts(Gender) <- contr.sum(2) # set factor-effects contrasts (i.e., dummy code 1 for male and -1 for female)

# Condition #
Condition <- with(data,convert(Q7.1)) # convert condition item
Condition[Condition==1] <- 'H' # recode 1 (high) --> H
Condition[Condition==2] <- 'L' # recode 2 (low) --> L
Condition <- factor(Condition,levels=c('H','L')) # designate Condition as factor object
contrasts(Condition) <- contr.sum(2) # set factor-effects contrasts (i.e., dummy code 1 for high and -1 for low)

# Subjective Commitment: Manipulation Check #     
SubjCommit.raw <- with(data,data.frame(Q3.3_1,Q3.3_2,Q3.3_3,Q3.3_4,Q3.3_5,Q3.3_6,Q3.3_7)) # read raw SubjCommit items
SubjCommit.data <- apply(SubjCommit.raw,2,convert) # convert SubjCommit items
SubjCommit.data[,4] <- 8-SubjCommit.data[,4] # reverse code 4th SubjCommit item
SubjCommit <- apply(SubjCommit.data,1,mean) # average SubjCommit across all items

# Exit Forgiveness #                 
Exit.raw <- with(data,data.frame(Q2.2_1,Q2.3_3,Q2.4_1,Q2.5_2,Q2.6_2,Q2.7_4,Q2.8_2,Q2.9_4,Q2.10_1,Q2.11_4,Q2.12_3,Q2.13_4)) # read Exit items
Exit.data <- apply(Exit.raw,2,convert) # convert Exit items
Exit <- apply(Exit.data,1,mean) # average Exit across all items

# Neglect Forgiveness #
Neglect.raw <- with(data,data.frame(Q2.2_4,Q2.3_2,Q2.4_2,Q2.5_4,Q2.6_1,Q2.7_2,Q2.8_4,Q2.9_2,Q2.10_4,Q2.11_1,Q2.12_1,Q2.13_1)) # read Neglect items
Neglect.data <- apply(Neglect.raw,2,convert) # convert Neglect items
Neglect <- apply(Neglect.data,1,mean) # average Neglect across all items

# Voice Forgiveness #
Voice.raw <- with(data,data.frame(Q2.2_2,Q2.3_10,Q2.4_3,Q2.5_1,Q2.6_4,Q2.7_1,Q2.8_3,Q2.9_3,Q2.10_2,Q2.11_2,Q2.12_4,Q2.13_3)) # read Voice items
Voice.data <- apply(Voice.raw,2,convert) # convert Voice items
Voice <- apply(Voice.data,1,mean) # average Voice across all items

# Loyalty Forgiveness #
Loyalty.raw <- with(data,data.frame(Q2.2_3,Q2.3_1,Q2.4_4,Q2.5_3,Q2.6_3,Q2.7_3,Q2.8_1,Q2.9_1,Q2.10_3,Q2.11_3, Q2.12_2,Q2.13_2)) #read Loyalty items
Loyalty.data <- apply(Loyalty.raw,2,convert) # convert Loyalty items
Loyalty <- apply(Loyalty.data,1,mean) # average Loyalty across all items

# Overall Forgiveness #
Overall <- apply(cbind(8-Exit,8-Neglect,Voice,Loyalty),1,mean) # calculate Overall by averaging the means of reversed Exit, reversed Neglect, Voice, and Loyalty

# Self Deception (SD) #
SD.raw <- data[,72:91] # read self deception items
SD.data <- apply(SD.raw,2,convert) # convert self deception items
SD.data[,seq(2,20,2)] <- 8-SD.data[,seq(2,20,2)] # reverse score even items
SD.data <- ifelse(SD.data>5,1,0) # extreme responses (6,7) = 1, else = 0
SD <- rowSums(SD.data) # number of extreme responses per subject
SD.c <- SD-mean(SD,na.rm=TRUE) # mean-center self deception (omitting missing scores in mean calculation)

# Impression Management (IM) #
IM.raw <- data[,92:111] # read impression management items
IM.data <- apply(IM.raw,2,convert) # convert impression management items
IM.data[,seq(1,19,2)] <- 8-IM.data[,seq(1,19,2)] # reverse score odd items
IM.data <- ifelse(IM.data>5,1,0) # extreme responses (6,7) = 1, else = 0
IM <- rowSums(IM.data) # number of extreme responses per subject
IM.c <- IM-mean(IM,na.rm=TRUE) # mean-center impression management (omitting missing scores in mean calculation)

# Formatted Data Table #
data.f <- data.frame(Gender,Condition,SubjCommit,Exit,Neglect,Voice,Loyalty,Overall,SD,IM) # fully formatted data table

#==========#
# Analyses #
#==========#

#------------------#
# Custom Functions #
#------------------#

# Source Table (Type III) #
aov3 <- function(model){ # model is an lm object
  tab <- drop1(model,~.,test='F')[-1,-3] # Type III SS including F-tests (delete <none> row and RSS column)
  tab[,3] <- tab[,2]/tab[,1] # replace AIC column with mean squares
  colnames(tab)[2:3] <- c('Sum Sq','Mean Sq') # rename SS and MS columns
  resid <- anova(model)[nrow(anova(model)),] # Df, SS, and MS for residuals (error)
  rbind(tab,resid) # add row for residuals
}

# Effect Size Analysis #
means <- function(model,alpha=.05,weighted=TRUE){ # model is an lm object, alpha=.05 by default, marginal means for each Condition weighted by gender proportions by default
  l <- length(model$coef) # number of parameters including intercept
  Y <- names(model$model)[1]
  if(l==2){ # t-test if only two parameters (intercept and Condition)
    results <- summary(model)$coef # regression estimates
    C <- matrix(c(results[1,1]+results[2,1],results[1,1]-results[2,1]),1,2) # means for high and low conditions
    dimnames(C) <- list('',c('H','L')) # row and column labels for C
    Diff <- C[1]-C[2] # difference in means (H-L)
    SE <- 2*results[2,2] # standard error of difference in means
    Df <- model$df # degrees of freedom
    ST <- cbind(Diff,SE,Df,results[2,3],results[2,4]) # t-test output
    dimnames(ST) <- list('',c('Diff (H-L)','SE','Df',colnames(results)[3:4])) # row and column labels for t-test output
    ES <- matrix(c(Diff,Diff+c(-1,1)*qt(1-alpha/2,Df)*SE),1,3) # mean difference and confidence interval output
    dimnames(ES) <- list('',c('Diff (H-L)','Lower','Upper')) # row and column labels for confidence interval output
    cat(paste('\n~~~~~~~~~~~~~~~~~~~~~~',paste(rep('~',nchar(Y)),sep='',collapse=''),sep=''))
    cat('\nTwo-Sample t-test for',paste(Y))
    cat(paste('\n~~~~~~~~~~~~~~~~~~~~~~',paste(rep('~',nchar(Y)),sep='',collapse=''),'\n',sep=''))
    cat('\nSample Means')
    cat('\n------------------\n')
    print(C)
    cat('\nt-test')
    cat('\n------------------------------------------------\n')
    print(ST)
    cat(paste('\n',(1-alpha)*100,'%',sep=''),'Confidence Interval')
    cat('\n---------------------------------\n')
    print(ES)
    cat('\n===========================================================\n')
    out <- list(C=C,ST=ST,ES=ES,meta=data.frame(H.Mean=C[1],L.Mean=C[2],Effect.Size=Diff,Standard.Error=SE)) # list of all outputs
  }
  else{ # ANOVA
    f <- which(attr(model$terms,'dataClasses')=='factor') # column positions of factors
    p <- c(1,f,l) # column positions of intercept, factors, and interaction
    X <- matrix(c(1,1,1,1,1,-1,1,-1,1,1,-1,-1,1,-1,-1,1),4,4) # unique design matrix without covariates (if any)
    B <- as.matrix(model$coef[p]) # parameter estimates without covariates (if any)
    C <- matrix(X%*%B,2,2) # (adjusted) cell means
    dimnames(C) <- model$xlevels # label rows and columns
    S <- vcov(model)[p,p] # covariance matrix of parameter estimates
    df <- model$df # error degrees of freedom
    n <- table(model$model[,f]) # cell sample sizes
    if(weighted==TRUE) { # for weighted marginal means of Condition by gender proportions
      w11 <- n[1,1]/(n[1,1]+n[2,1]) # weighting constant for cell mean of male / high
      w21 <- n[2,1]/(n[1,1]+n[2,1]) # weighting constant for cell mean of female / high
      w12 <- n[1,2]/(n[1,2]+n[2,2]) # weighting constant for cell mean of male / low
      w22 <- n[2,2]/(n[1,2]+n[2,2]) # weighting constant for cell mean of female / low
      wa <- (n[1,1]-n[2,1])/(n[1,1]+n[2,1])-(n[1,2]-n[2,2])/(n[1,2]+n[2,2]) # weighting constant for a
      wab <- (n[1,1]-n[2,1])/(n[1,1]+n[2,1])+(n[1,2]-n[2,2])/(n[1,2]+n[2,2]) # weighting constant for ab
      W <- 'Weighted'
    }
    if(weighted==FALSE) { # for unweighted marginal means of Condition
      w11 <- w21 <- w12 <- w22 <- .5 # equal weights
      wa <- wab <- 0 # zero weights
      W <- 'Unweighted'
    }
    MH <- w11*C[1,1]+w21*C[2,1] # marginal mean for high
    ML <- w12*C[1,2]+w22*C[2,2] # marginal mean for low
    M <- matrix(c(MH,ML),1,2) # marginal means
    dimnames(M) <- list(W,c('H','L')) # condition labels
    HM <- c(MH,C[1,1],C[2,1]) # means for high
    LM <- c(ML,C[1,2],C[2,2]) # means for low
    DV_D <- 2*B[3]+wa*B[2]+wab*B[4] # difference in marginal means (H-L)
    DV_SE <- sqrt(4*S[3,3]+wa^2*S[2,2]+wab^2*S[4,4]+4*wa*S[3,2]+4*wab*S[3,4]+2*wa*wab*S[2,4]) # standard error of difference in marginal means
    CI_D <- DV_D + c(-1,1)*qt(1-alpha/6,df)*DV_SE # Bonferroni confidence interval for difference in marginal means
    DV_D_M <- 2*(B[3]+B[4]) # difference in means for males
    DV_SE_M <- 2*sqrt(S[3,3]+S[4,4]+2*S[3,4]) # standard error of difference in means for males
    CI_D_M <- DV_D_M + c(-1,1)*qt(1-alpha/6,df)*DV_SE_M # Bonferroni confidence interval for difference in means for males
    DV_D_F <- 2*(B[3]-B[4]) # difference in means for females
    DV_SE_F <- 2*sqrt(S[3,3]+S[4,4]-2*S[3,4]) # standard error of difference in means for females
    CI_D_F <- DV_D_F + c(-1,1)*qt(1-alpha/6,df)*DV_SE_F # Bonferroni confidence interval for difference in means for females
    CI <- rbind(CI_D,CI_D_M,CI_D_F) # output of all confidence intervals
    DV <- c(DV_D,DV_D_M,DV_D_F) # output of all differences in means
    SE <- c(DV_SE,DV_SE_M,DV_SE_F) # output of all standard errors
    ES <- cbind(DV,CI) # output of all differences in means and standard errors
    dimnames(ES) <- list(c(paste(W,'Marginal'),'Male','Female'),c('Diff (H-L)',' Lower','Upper')) # row and column labels of ES output
    ST <- aov3(model) # Type III ANOVA output
    cat(paste('\n~~~~~~~~~~~~~~~',paste(rep('~',nchar(Y)),sep='',collapse=''),sep=''))
    if(l==4) cat('\n2x2 ANOVA for',paste(Y)) 
    else cat('\n2x2 ANCOVA for',paste(Y))
    cat(paste('\n~~~~~~~~~~~~~~~',paste(rep('~',nchar(Y)),sep='',collapse=''),'\n',sep=''))
    cat('\nSource Table (Type III)')
    cat('\n-------------------------------------------------------\n')
    print(ST,signif.legend=FALSE,signif.stars=FALSE)
    cat('\n')
    if(l==4) cat('Cell Means')
    else cat('Adjusted Cell Means')
    cat('\n------------------------\n')
    print(C)
    cat('\n')
    if(l==4) cat('Marginal Means')
    else cat('Adjusted Marginal Means')
    cat('\n----------------------------\n')
    print(M)
    cat(paste('\n',(1-alpha)*100,'%',sep=''),'Bonferroni Confidence Intervals')
    cat('\n--------------------------------------------------\n')
    print(ES)
    cat('\n===========================================================\n')
    out <- list(ST=ST,C=C,M=M,ES=ES,meta=data.frame(H.Mean=HM,L.Mean=LM,Effect.Size=DV,Standard.Error=SE)) # list of all outputs
  }
}

# Mediation Analysis (Bias-Corrected Bootstrap Method) #
mediation.bs <- function(data,Y,M,X,r,alpha=.05){ # Y = dependent variable, M = mediator, X = independent variable, r = number of replications, alpha=.05 by default
  a <- lm(data[,M]~data[,X])$coef[2] # original a in M = aX
  b <- lm(data[,Y]~data[,M]+data[,X])$coef[2] # original b in Y = bM + cX
  ab0 <- a*b # a*b (estimated mediation effect based on original sample)
  ab <- rep(0,r) # allocate a*b vector for r replications
  n <- nrow(data) # sample size
  for(i in 1:r){ # bootstrap replications
    data.bs <- data[sample(n,n,replace=TRUE),] # draw bootstrap sample
    a <- lm(data.bs[,M]~data.bs[,X])$coef[2] # replicate a
    b <- lm(data.bs[,Y]~data.bs[,M]+data.bs[,X])$coef[2] # replicate b
    ab[i] <- a*b # replicate a*b
  }
  mean.ab <- mean(ab) # bootstrap estimate of a*b (mean of all replicated a*b)
  sd.ab <- sd(ab) # bootstrap estimate of standard error of a*b (standard deviation of all replicated a*b)
  bc <- 2*qnorm(mean(ab<ab0)) # measure of bias in terms of z quantile (away from center of 0)
  zl <- bc+qnorm(alpha/2) # lower quantile = bias + z(.05)
  zu <- bc+qnorm(1-alpha/2) # upper quantile = bias + z(.95)
  CI <- t(matrix(quantile(ab,c(pnorm(zl),pnorm(zu))))) # bias-corrected confidence interval
  results <- cbind(mean.ab,CI) # output of estimated a*b and confidence interval
  dimnames(results) <- list('a*b',c('Estimate','Lower','Upper')) # row and column labels for output
  cat(paste('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~',paste(rep('~',nchar(Y)),sep='',collapse=''),'\n',sep=''))
  cat(paste('Mediation of Condition on',paste(Y),'by SubjCommit'))
  cat(paste('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~',paste(rep('~',nchar(Y)),sep='',collapse=''),'\n',sep=''))
  cat(paste('\n',(1-alpha)*100,'%',sep=''),'Bias-Corrected Confidence Interval',paste('(',r,' replications)',sep=''))
  cat('\n-----------------------------------------------------------\n')
  print(results)
  cat('\n===========================================================\n')
  out <- list(results=results,meta=data.frame(H.Mean=NA,L.Mean=NA,Effect.Size=mean.ab,Standard.Error=sd.ab)) # list of all outputs
}

#---------#
# Results #
#---------#

# Subjective Commitment (Manipulation Check) #
SubjCommit.t <- means(lm(SubjCommit~Condition))
SubjCommit.ANOVA <- means(lm(SubjCommit~Gender*Condition))
SubjCommit.ANCOVA <- means(lm(SubjCommit~SD.c+IM.c+Gender*Condition))
SubjCommit.meta <- rbind(SubjCommit.t$meta,SubjCommit.ANOVA$meta,SubjCommit.ANCOVA$meta)
rownames(SubjCommit.meta) <- c('High-Low (t)','High-Low (ANOVA)','High-Low: Male (ANOVA)','High-Low: Female (ANOVA)','High-Low (ANCOVA)','High-Low: Male (ANCOVA)','High-Low: Female (ANCOVA)')

# Exit Forgiveness #
Exit.t <- means(lm(Exit~Condition))
Exit.ANOVA <- means(lm(Exit~Gender*Condition))
Exit.ANCOVA <- means(lm(Exit~SD.c+IM.c+Gender*Condition))
Exit.med <- mediation.bs(data.f,Y='Exit',M='SubjCommit',X='Condition',r)
Exit.meta <- rbind(Exit.t$meta,Exit.ANOVA$meta,Exit.ANCOVA$meta,Exit.med$meta)
rownames(Exit.meta) <- c('High-Low (t)','High-Low (ANOVA)','High-Low: Male (ANOVA)','High-Low: Female (ANOVA)','High-Low (ANCOVA)','High-Low: Male (ANCOVA)','High-Low: Female (ANCOVA)','Mediation (a*b)')

# Neglect Forgiveness #
Neglect.t <- means(lm(Neglect~Condition))
Neglect.ANOVA <- means(lm(Neglect~Gender*Condition))
Neglect.ANCOVA <- means(lm(Neglect~SD.c+IM.c+Gender*Condition))
Neglect.med <- mediation.bs(data.f,Y='Neglect',M='SubjCommit',X='Condition',r)
Neglect.meta <- rbind(Neglect.t$meta,Neglect.ANOVA$meta,Neglect.ANCOVA$meta,Neglect.med$meta)
rownames(Neglect.meta) <- c('High-Low (t)','High-Low (ANOVA)','High-Low: Male (ANOVA)','High-Low: Female (ANOVA)','High-Low (ANCOVA)','High-Low: Male (ANCOVA)','High-Low: Female (ANCOVA)','Mediation (a*b)')

# Voice Forgiveness #
Voice.t <- means(lm(Voice~Condition))
Voice.ANOVA <- means(lm(Voice~Gender*Condition))
Voice.ANCOVA <- means(lm(Voice~SD.c+IM.c+Gender*Condition))
Voice.meta <- rbind(Voice.t$meta,Voice.ANOVA$meta,Voice.ANCOVA$meta)
rownames(Voice.meta) <- c('High-Low (t)','High-Low (ANOVA)','High-Low: Male (ANOVA)','High-Low: Female (ANOVA)','High-Low (ANCOVA)','High-Low: Male (ANCOVA)','High-Low: Female (ANCOVA)')

# Loyalty Forgiveness #
Loyalty.t <- means(lm(Loyalty~Condition))
Loyalty.ANOVA <- means(lm(Loyalty~Gender*Condition))
Loyalty.ANCOVA <- means(lm(Loyalty~SD.c+IM.c+Gender*Condition))
Loyalty.meta <- rbind(Loyalty.t$meta,Loyalty.ANOVA$meta,Loyalty.ANCOVA$meta)
rownames(Loyalty.meta) <- c('High-Low (t)','High-Low (ANOVA)','High-Low: Male (ANOVA)','High-Low: Female (ANOVA)','High-Low (ANCOVA)','High-Low: Male (ANCOVA)','High-Low: Female (ANCOVA)')

# Overall Forgiveness #
Overall.t <- means(lm(Overall~Condition))
Overall.ANOVA <- means(lm(Overall~Gender*Condition))
Overall.ANCOVA <- means(lm(Overall~SD.c+IM.c+Gender*Condition))
Overall.meta <- rbind(Overall.t$meta,Overall.ANOVA$meta,Overall.ANCOVA$meta)
rownames(Overall.meta) <- c('High-Low (t)','High-Low (ANOVA)','High-Low: Male (ANOVA)','High-Low: Female (ANOVA)','High-Low (ANCOVA)','High-Low: Male (ANCOVA)','High-Low: Female (ANCOVA)')

meta <- list(SubjCommit=SubjCommit.meta,Exit=Exit.meta,Neglect=Neglect.meta,Voice=Voice.meta,Loyalty=Loyalty.meta,Overall=Overall.meta)

}

file.csv <- list.files(pattern="*.csv")[1] # name of csv file in source folder
  # if there are multiple csv files in folder, use first one according to alphabetical order
author <- strsplit(file.csv,split='_',fixed=TRUE)[[1]][1] # extract author's name
data <- read.csv(file.csv[1],header=FALSE) # read csv file

data2 <- data # duplicate data
colnames(data2) <- apply(data2[1,],2,as.character) # label columns
data2 <- data2[-(1:2),] # delete first two irrelevant rows
Exclude <- convert(data2$Exclude) # exclusion criteria
if(max(Exclude)==2) m <- 2 else m <- 1 # if no subjects are Exclude=2, run analyses just once
rm(Exclude,data2) # remove Exclude and data2 from environment

for(i in 1:m){

if(i==1) include <- 0 # first run with only 0
if(i==2) include <- c(0,2) # second run with 0 and 2

file.out <- paste(author,'_output_',paste(include,collapse=''),'.txt',sep='') # start writing named output file
sink(file.out) # start writing all console output to named text file
cat('===========================\n')
cat('FINKEL REPLICATION ANALYSES\n')
cat('===========================\n')
cat('\nInclusion\n')
cat('---------\n')
cat('0 = Include Subject\n')
cat('1 = Exclude Subject\n')
cat('2 = Maybe Include Subject\n')
cat('\nEffect Size\n')
cat('-----------\n')
cat('Difference in Means between High and Low Priming Conditions\n')
cat('\nDependent Variables\n')
cat('-------------------\n')
cat('Subjective Commitment (SubjCommit)\n') 
cat('Exit Forgiveness (Exit)\n')
cat('Neglect Forgiveness (Neglect)\n')
cat('Voice Forgiveness (Voice)\n')
cat('Loyalty Forgiveness (Loyalty)\n')
cat('Overall Forgiveness (Overall)\n')
cat('\nIndependent Variable\n')
cat('--------------------\n')
cat('Priming Condition (Condition): High (H), Low (L)\n')
cat('\nModerator\n')
cat('---------\n')
cat('Gender (Gender): Male (M), Female (F)\n')
cat('\nCovariates\n') 
cat('----------\n')
cat('Self Deception (SD)\n')
cat('Impression Management (IM)\n')
cat('\nAnalyses Output for Each Dependent Variable\n')
cat('-------------------------------------------\n')
cat('Two-sample t-test\n')
cat(' (1) Sample Means for High and Low\n')
cat(' (2) t-test of Difference in Means\n')
cat(' (3) 95% Confidence Intervals for Difference in Means\n')
cat('\n2x2 ANOVA (Gender x Condition) \n')
cat(' (1) Source Table (Type III Sums of Squares) \n')
cat(' (2) Cell means (Male/High, Male/Low, Female/High, Female/Low) \n')
cat(' (3) Weighted Marginal Means across Gender for High and Low\n')
cat('     *marginal means weighted by Gender proportions in each Condition\n')
cat('     *equal to sample means for two-sample t-test\n')
cat(' (4) 95% Bonferroni Confidence Intervals (family of 3) for Difference in: \n')
cat('      (i)   Weighted Marginal Means across Gender for High and Low\n')
cat('      (ii)  Male Means for High and Low\n')
cat('      (iii) Female Means for High and Low\n')
cat('\n2x2 ANCOVA (Gender x Condition with SD and IM covariates) \n')
cat(' (1) Source Table (Type III Sums of Squares) \n')
cat(' (2) Adjusted Cell means (Male/High, Male/Low, Female/High, Female/Low) \n')
cat('     *adjusted mean = mean of dependent variable at means of covariates \n')
cat(' (3) Weighted Marginal Adjusted Means across Gender for High and Low\n')
cat('     *marginal adjusted means weighted by Gender proportions in each Condition\n')
cat(' (4) 95% Bonferroni Confidence Intervals (family of 3) for Difference in: \n')
cat('      (i)   Weighted Marginal Adjusted Means across Gender for High and Low\n')
cat('      (ii)  Male Adjusted Means for High and Low\n')
cat('      (iii) Female Adjusted Means for High and Low\n')
cat('\nMediation Analysis (for Exit and Neglect only)\n')
cat(" (1) 95% Bias-Corrected Confidence Interval for Subjective Commitment's \n")
cat('     Mediation of the Effect of Condition on Exit or Neglect Forgiveness\n')
cat('     *constructed using bootstrap with 10,000 replications\n')
cat('\n*******************************************************************************\n')
cat(paste('\n============',paste(rep('=',nchar(author)),sep='',collapse=''),'\n',sep=''))
cat(paste('RESULTS FOR',toupper(author)))
cat(paste('\n============',paste(rep('=',nchar(author)),sep='',collapse=''),'\n',sep=''))
cat(paste('\n*Inclusion:',paste(include,collapse=' & '),'\n'))
Finkel.Analyze(data,include) # run function with data
sink() # end writing output

}

rm(author,include,file.csv,file.out,i,m,Finkel.Analyze,data,convert) # clean environment

cat('\nDone! Find results in folder.')
