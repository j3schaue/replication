setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../raw/rrr_cheung")

library(dplyr)

#Convert Function
convert <- function(x) as.numeric(as.character(x)) # function that converts a column of Qualtrix output to useable form in R

#Analysis Function
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
    #out <- list(C=C,ST=ST,ES=ES,meta=data.frame(H.Mean=C[1],L.Mean=C[2],Effect.Size=Diff,Standard.Error=SE)) # list of all outputs
    out <- c(ST) # list of all outputs
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
include <- 0

#Read Data
ay <- read.csv("./Data/Aykutoglu_Data.csv", fileEncoding="latin1")
ay <- ay[-(1:2),] # delete first two irrelevant rows
ayExclude <- convert(ay$Exclude) # exclusion criteria
ay <- ay[which(ayExclude %in% include),] # include specified observations
aycond <- with(ay,convert(Q7.1)) # convert condition item
aycond[aycond==1] <- 'H' # recode 1 (high) --> H
aycond[aycond==2] <- 'L' # recode 2 (low) --> L
aycond <- factor(aycond,levels=c('H','L')) # designate Condition as factor object
contrasts(aycond) <- contr.sum(2) # set factor-effects contrasts (i.e., dummy code 1 for high and -1 for low)
ayExit.raw <- with(ay,data.frame(Q2.2_1,Q2.3_3,Q2.4_1,Q2.5_2,Q2.6_2,Q2.7_4,Q2.8_2,Q2.9_4,Q2.10_1,Q2.11_4,Q2.12_3,Q2.13_4)) # read Exit items
ayExit.data <- apply(ayExit.raw,2,convert) # convert Exit items
ayExit <- apply(ayExit.data,1,mean) # average Exit across all items
ayNeglect.raw <- with(ay,data.frame(Q2.2_4,Q2.3_2,Q2.4_2,Q2.5_4,Q2.6_1,Q2.7_2,Q2.8_4,Q2.9_2,Q2.10_4,Q2.11_1,Q2.12_1,Q2.13_1)) # read Neglect items
ayNeglect.data <- apply(ayNeglect.raw,2,convert) # convert Neglect items
ayNeglect <- apply(ayNeglect.data,1,mean) # average Neglect across all items
ayVoice.raw <- with(ay,data.frame(Q2.2_2,Q2.3_10,Q2.4_3,Q2.5_1,Q2.6_4,Q2.7_1,Q2.8_3,Q2.9_3,Q2.10_2,Q2.11_2,Q2.12_4,Q2.13_3)) # read Voice items
ayVoice.data <- apply(ayVoice.raw,2,convert) # convert Voice items
ayVoice <- apply(ayVoice.data,1,mean) # average Voice across all items
ayLoyalty.raw <- with(ay,data.frame(Q2.2_3,Q2.3_1,Q2.4_4,Q2.5_3,Q2.6_3,Q2.7_3,Q2.8_1,Q2.9_1,Q2.10_3,Q2.11_3, Q2.12_2,Q2.13_2)) #read Loyalty items
ayLoyalty.data <- apply(ayLoyalty.raw,2,convert) # convert Loyalty items
ayLoyalty <- apply(ayLoyalty.data,1,mean) # average Loyalty across all items

bredow <- read.csv("./Data/Bredow_data.csv", fileEncoding="latin1")
bredow <- bredow[-(1:2),]
bExclude <- convert(bredow$Exclude)
bredow <- bredow[which(bExclude %in% include),]
bcond <- with(bredow,convert(Q7.1))
bcond[bcond==1] <- 'H'
bcond[bcond==2] <- 'L'
bcond <- factor(bcond,levels=c('H','L'))
contrasts(bcond) <- contr.sum(2)
bExit.raw <- with(bredow,data.frame(Q2.2_1,Q2.3_3,Q2.4_1,Q2.5_2,Q2.6_2,Q2.7_4,Q2.8_2,Q2.9_4,Q2.10_1,Q2.11_4,Q2.12_3,Q2.13_4)) # read Exit items
bExit.data <- apply(bExit.raw,2,convert)
bExit <- apply(bExit.data,1,mean)
bNeglect.raw <- with(bredow,data.frame(Q2.2_4,Q2.3_2,Q2.4_2,Q2.5_4,Q2.6_1,Q2.7_2,Q2.8_4,Q2.9_2,Q2.10_4,Q2.11_1,Q2.12_1,Q2.13_1)) # read Neglect items
bNeglect.data <- apply(bNeglect.raw,2,convert)
bNeglect <- apply(bNeglect.data,1,mean)
bVoice.raw <- with(bredow,data.frame(Q2.2_2,Q2.3_10,Q2.4_3,Q2.5_1,Q2.6_4,Q2.7_1,Q2.8_3,Q2.9_3,Q2.10_2,Q2.11_2,Q2.12_4,Q2.13_3)) # read Voice items
bVoice.data <- apply(bVoice.raw,2,convert)
bVoice <- apply(bVoice.data,1,mean)
bLoyalty.raw <- with(bredow,data.frame(Q2.2_3,Q2.3_1,Q2.4_4,Q2.5_3,Q2.6_3,Q2.7_3,Q2.8_1,Q2.9_1,Q2.10_3,Q2.11_3, Q2.12_2,Q2.13_2)) #read Loyalty items
bLoyalty.data <- apply(bLoyalty.raw,2,convert)
bLoyalty <- apply(bLoyalty.data,1,mean)

cap <- read.csv("./Data/Caprariello_data.csv", fileEncoding="latin1")
cap <- cap[-(1:2),]
capExclude <- convert(cap$Exclude)
cap <- cap[which(capExclude %in% include),]
capcond <- with(cap,convert(Q7.1))
capcond[capcond==1] <- 'H'
capcond[capcond==2] <- 'L'
capcond <- factor(capcond,levels=c('H','L'))
contrasts(capcond) <- contr.sum(2)
capExit.raw <- with(cap,data.frame(Q2.2_1,Q2.3_3,Q2.4_1,Q2.5_2,Q2.6_2,Q2.7_4,Q2.8_2,Q2.9_4,Q2.10_1,Q2.11_4,Q2.12_3,Q2.13_4)) # read Exit items
capExit.data <- apply(capExit.raw,2,convert)
capExit <- apply(capExit.data,1,mean)
capNeglect.raw <- with(cap,data.frame(Q2.2_4,Q2.3_2,Q2.4_2,Q2.5_4,Q2.6_1,Q2.7_2,Q2.8_4,Q2.9_2,Q2.10_4,Q2.11_1,Q2.12_1,Q2.13_1)) # read Neglect items
capNeglect.data <- apply(capNeglect.raw,2,convert)
capNeglect <- apply(capNeglect.data,1,mean)
capVoice.raw <- with(cap,data.frame(Q2.2_2,Q2.3_10,Q2.4_3,Q2.5_1,Q2.6_4,Q2.7_1,Q2.8_3,Q2.9_3,Q2.10_2,Q2.11_2,Q2.12_4,Q2.13_3)) # read Voice items
capVoice.data <- apply(capVoice.raw,2,convert)
capVoice <- apply(capVoice.data,1,mean)
capLoyalty.raw <- with(cap,data.frame(Q2.2_3,Q2.3_1,Q2.4_4,Q2.5_3,Q2.6_3,Q2.7_3,Q2.8_1,Q2.9_1,Q2.10_3,Q2.11_3, Q2.12_2,Q2.13_2)) #read Loyalty items
capLoyalty.data <- apply(capLoyalty.raw,2,convert)
capLoyalty <- apply(capLoyalty.data,1,mean)

carcedo <- read.csv("./Data/Carcedo_data.csv", fileEncoding="latin1")
carcedo <- carcedo[-(1:2),]
carExclude <- convert(carcedo$Exclude)
carcedo <- carcedo[which(carExclude %in% include),]
carcond <- with(carcedo,convert(Q7.1))
carcond[carcond==1] <- 'H'
carcond[carcond==2] <- 'L'
carcond <- factor(carcond,levels=c('H','L'))
contrasts(carcond) <- contr.sum(2)
carExit.raw <- with(carcedo,data.frame(Q2.2_1,Q2.3_3,Q2.4_1,Q2.5_2,Q2.6_2,Q2.7_4,Q2.8_2,Q2.9_4,Q2.10_1,Q2.11_4,Q2.12_3,Q2.13_4)) # read Exit items
carExit.data <- apply(carExit.raw,2,convert)
carExit <- apply(carExit.data,1,mean)
carNeglect.raw <- with(carcedo,data.frame(Q2.2_4,Q2.3_2,Q2.4_2,Q2.5_4,Q2.6_1,Q2.7_2,Q2.8_4,Q2.9_2,Q2.10_4,Q2.11_1,Q2.12_1,Q2.13_1)) # read Neglect items
carNeglect.data <- apply(carNeglect.raw,2,convert)
carNeglect <- apply(carNeglect.data,1,mean)
carVoice.raw <- with(carcedo,data.frame(Q2.2_2,Q2.3_10,Q2.4_3,Q2.5_1,Q2.6_4,Q2.7_1,Q2.8_3,Q2.9_3,Q2.10_2,Q2.11_2,Q2.12_4,Q2.13_3)) # read Voice items
carVoice.data <- apply(carVoice.raw,2,convert)
carVoice <- apply(carVoice.data,1,mean)
carLoyalty.raw <- with(carcedo,data.frame(Q2.2_3,Q2.3_1,Q2.4_4,Q2.5_3,Q2.6_3,Q2.7_3,Q2.8_1,Q2.9_1,Q2.10_3,Q2.11_3, Q2.12_2,Q2.13_2)) #read Loyalty items
carLoyalty.data <- apply(carLoyalty.raw,2,convert)
carLoyalty <- apply(carLoyalty.data,1,mean)

carson <- read.csv("./Data/Carson_data.csv", fileEncoding="latin1")
carson <- carson[-(1:2),]
carExclude <- convert(carson$Exclude)
carson <- carson[which(carExclude %in% include),]
carscond <- with(carson,convert(Q7.1))
carscond[carscond==1] <- 'H'
carscond[carscond==2] <- 'L'
carscond <- factor(carscond,levels=c('H','L'))
contrasts(carscond) <- contr.sum(2)
carsExit.raw <- with(carson,data.frame(Q2.2_1,Q2.3_3,Q2.4_1,Q2.5_2,Q2.6_2,Q2.7_4,Q2.8_2,Q2.9_4,Q2.10_1,Q2.11_4,Q2.12_3,Q2.13_4)) # read Exit items
carsExit.data <- apply(carsExit.raw,2,convert)
carsExit <- apply(carsExit.data,1,mean)
carsNeglect.raw <- with(carson,data.frame(Q2.2_4,Q2.3_2,Q2.4_2,Q2.5_4,Q2.6_1,Q2.7_2,Q2.8_4,Q2.9_2,Q2.10_4,Q2.11_1,Q2.12_1,Q2.13_1)) # read Neglect items
carsNeglect.data <- apply(carsNeglect.raw,2,convert)
carsNeglect <- apply(carsNeglect.data,1,mean)
carsVoice.raw <- with(carson,data.frame(Q2.2_2,Q2.3_10,Q2.4_3,Q2.5_1,Q2.6_4,Q2.7_1,Q2.8_3,Q2.9_3,Q2.10_2,Q2.11_2,Q2.12_4,Q2.13_3)) # read Voice items
carsVoice.data <- apply(carsVoice.raw,2,convert)
carsVoice <- apply(carsVoice.data,1,mean)
carsLoyalty.raw <- with(carson,data.frame(Q2.2_3,Q2.3_1,Q2.4_4,Q2.5_3,Q2.6_3,Q2.7_3,Q2.8_1,Q2.9_1,Q2.10_3,Q2.11_3, Q2.12_2,Q2.13_2)) #read Loyalty items
carsLoyalty.data <- apply(carsLoyalty.raw,2,convert)
carsLoyalty <- apply(carsLoyalty.data,1,mean)

cheung <- read.csv("./Data/Cheung_data.csv", fileEncoding="latin1")
cheung <- cheung[-(1:2),] 
chExclude <- convert(cheung$Exclude) 
cheung <- cheung[which(chExclude %in% include),]
chcond <- with(cheung,convert(Q7.1))
chcond[chcond==1] <- 'H'
chcond[chcond==2] <- 'L'
chcond <- factor(chcond,levels=c('H','L'))
contrasts(chcond) <- contr.sum(2)
chExit.raw <- with(cheung,data.frame(Q2.2_1,Q2.3_3,Q2.4_1,Q2.5_2,Q2.6_2,Q2.7_4,Q2.8_2,Q2.9_4,Q2.10_1,Q2.11_4,Q2.12_3,Q2.13_4)) # read Exit items
chExit.data <- apply(chExit.raw,2,convert)
chExit <- apply(chExit.data,1,mean)
chNeglect.raw <- with(cheung,data.frame(Q2.2_4,Q2.3_2,Q2.4_2,Q2.5_4,Q2.6_1,Q2.7_2,Q2.8_4,Q2.9_2,Q2.10_4,Q2.11_1,Q2.12_1,Q2.13_1)) # read Neglect items
chNeglect.data <- apply(chNeglect.raw,2,convert)
chNeglect <- apply(chNeglect.data,1,mean)
chVoice.raw <- with(cheung,data.frame(Q2.2_2,Q2.3_10,Q2.4_3,Q2.5_1,Q2.6_4,Q2.7_1,Q2.8_3,Q2.9_3,Q2.10_2,Q2.11_2,Q2.12_4,Q2.13_3)) # read Voice items
chVoice.data <- apply(chVoice.raw,2,convert)
chVoice <- apply(chVoice.data,1,mean)
chLoyalty.raw <- with(cheung,data.frame(Q2.2_3,Q2.3_1,Q2.4_4,Q2.5_3,Q2.6_3,Q2.7_3,Q2.8_1,Q2.9_1,Q2.10_3,Q2.11_3, Q2.12_2,Q2.13_2)) #read Loyalty items
chLoyalty.data <- apply(chLoyalty.raw,2,convert)
chLoyalty <- apply(chLoyalty.data,1,mean)

cobb <- read.csv("./Data/Cobb_data.csv", fileEncoding="latin1")
cobb <- cobb[-(1:2),]
coExclude <- convert(cobb$Exclude)
cobb <- cobb[which(coExclude %in% include),]
cocond <- with(cobb,convert(Q7.1))
cocond[cocond==1] <- 'H'
cocond[cocond==2] <- 'L'
cocond <- factor(cocond,levels=c('H','L'))
contrasts(cocond) <- contr.sum(2)
coExit.raw <- with(cobb,data.frame(Q2.2_1,Q2.3_3,Q2.4_1,Q2.5_2,Q2.6_2,Q2.7_4,Q2.8_2,Q2.9_4,Q2.10_1,Q2.11_4,Q2.12_3,Q2.13_4)) # read Exit items
coExit.data <- apply(coExit.raw,2,convert)
coExit <- apply(coExit.data,1,mean)
coNeglect.raw <- with(cobb,data.frame(Q2.2_4,Q2.3_2,Q2.4_2,Q2.5_4,Q2.6_1,Q2.7_2,Q2.8_4,Q2.9_2,Q2.10_4,Q2.11_1,Q2.12_1,Q2.13_1)) # read Neglect items
coNeglect.data <- apply(coNeglect.raw,2,convert)
coNeglect <- apply(coNeglect.data,1,mean)
coVoice.raw <- with(cobb,data.frame(Q2.2_2,Q2.3_10,Q2.4_3,Q2.5_1,Q2.6_4,Q2.7_1,Q2.8_3,Q2.9_3,Q2.10_2,Q2.11_2,Q2.12_4,Q2.13_3)) # read Voice items
coVoice.data <- apply(coVoice.raw,2,convert)
coVoice <- apply(coVoice.data,1,mean)
coLoyalty.raw <- with(cobb,data.frame(Q2.2_3,Q2.3_1,Q2.4_4,Q2.5_3,Q2.6_3,Q2.7_3,Q2.8_1,Q2.9_1,Q2.10_3,Q2.11_3, Q2.12_2,Q2.13_2)) #read Loyalty items
coLoyalty.data <- apply(coLoyalty.raw,2,convert)
coLoyalty <- apply(coLoyalty.data,1,mean)

collin <- read.csv("./Data/Collins_Data.csv", fileEncoding="latin1")
collin <- collin[-(1:2),]
collinExclude <- convert(collin$Exclude)
collin <- collin[which(collinExclude %in% include),]
collincond <- with(collin,convert(Q7.1))
collincond[collincond==1] <- 'H'
collincond[collincond==2] <- 'L'
collincond <- factor(collincond,levels=c('H','L'))
contrasts(collincond) <- contr.sum(2)
collinExit.raw <- with(collin,data.frame(Q2.2_1,Q2.3_3,Q2.4_1,Q2.5_2,Q2.6_2,Q2.7_4,Q2.8_2,Q2.9_4,Q2.10_1,Q2.11_4,Q2.12_3,Q2.13_4)) # read Exit items
collinExit.data <- apply(collinExit.raw,2,convert)
collinExit <- apply(collinExit.data,1,mean)
collinNeglect.raw <- with(collin,data.frame(Q2.2_4,Q2.3_2,Q2.4_2,Q2.5_4,Q2.6_1,Q2.7_2,Q2.8_4,Q2.9_2,Q2.10_4,Q2.11_1,Q2.12_1,Q2.13_1)) # read Neglect items
collinNeglect.data <- apply(collinNeglect.raw,2,convert)
collinNeglect <- apply(collinNeglect.data,1,mean)
collinVoice.raw <- with(collin,data.frame(Q2.2_2,Q2.3_10,Q2.4_3,Q2.5_1,Q2.6_4,Q2.7_1,Q2.8_3,Q2.9_3,Q2.10_2,Q2.11_2,Q2.12_4,Q2.13_3)) # read Voice items
collinVoice.data <- apply(collinVoice.raw,2,convert)
collinVoice <- apply(collinVoice.data,1,mean)
collinLoyalty.raw <- with(collin,data.frame(Q2.2_3,Q2.3_1,Q2.4_4,Q2.5_3,Q2.6_3,Q2.7_3,Q2.8_1,Q2.9_1,Q2.10_3,Q2.11_3, Q2.12_2,Q2.13_2)) #read Loyalty items
collinLoyalty.data <- apply(collinLoyalty.raw,2,convert)
collinLoyalty <- apply(collinLoyalty.data,1,mean)

dido <- read.csv("./Data/Didonato_data.csv", fileEncoding="latin1")
dido <- dido[-(1:2),]
dExclude <- convert(dido$Exclude)
dido <- dido[which(dExclude %in% include),]
dcond <- with(dido,convert(Q7.1))
dcond[dcond==1] <- 'H'
dcond[dcond==2] <- 'L'
dcond <- factor(dcond,levels=c('H','L'))
contrasts(dcond) <- contr.sum(2)
dExit.raw <- with(dido,data.frame(Q2.2_1,Q2.3_3,Q2.4_1,Q2.5_2,Q2.6_2,Q2.7_4,Q2.8_2,Q2.9_4,Q2.10_1,Q2.11_4,Q2.12_3,Q2.13_4)) # read Exit items
dExit.data <- apply(dExit.raw,2,convert)
dExit <- apply(dExit.data,1,mean)
dNeglect.raw <- with(dido,data.frame(Q2.2_4,Q2.3_2,Q2.4_2,Q2.5_4,Q2.6_1,Q2.7_2,Q2.8_4,Q2.9_2,Q2.10_4,Q2.11_1,Q2.12_1,Q2.13_1)) # read Neglect items
dNeglect.data <- apply(dNeglect.raw,2,convert)
dNeglect <- apply(dNeglect.data,1,mean)
dVoice.raw <- with(dido,data.frame(Q2.2_2,Q2.3_10,Q2.4_3,Q2.5_1,Q2.6_4,Q2.7_1,Q2.8_3,Q2.9_3,Q2.10_2,Q2.11_2,Q2.12_4,Q2.13_3)) # read Voice items
dVoice.data <- apply(dVoice.raw,2,convert)
dVoice <- apply(dVoice.data,1,mean)
dLoyalty.raw <- with(dido,data.frame(Q2.2_3,Q2.3_1,Q2.4_4,Q2.5_3,Q2.6_3,Q2.7_3,Q2.8_1,Q2.9_1,Q2.10_3,Q2.11_3, Q2.12_2,Q2.13_2)) #read Loyalty items
dLoyalty.data <- apply(dLoyalty.raw,2,convert)
dLoyalty <- apply(dLoyalty.data,1,mean)

fug <- read.csv("./Data/Fuglestad_data.csv", fileEncoding="latin1")
fug <- fug[-(1:2),]
fExclude <- convert(fug$Exclude)
fug <- fug[which(fExclude %in% include),]
fcond <- with(fug,convert(Q7.1))
fcond[fcond==1] <- 'H'
fcond[fcond==2] <- 'L'
fcond <- factor(fcond,levels=c('H','L'))
contrasts(fcond) <- contr.sum(2)
fExit.raw <- with(fug,data.frame(Q2.2_1,Q2.3_3,Q2.4_1,Q2.5_2,Q2.6_2,Q2.7_4,Q2.8_2,Q2.9_4,Q2.10_1,Q2.11_4,Q2.12_3,Q2.13_4)) # read Exit items
fExit.data <- apply(fExit.raw,2,convert)
fExit <- apply(fExit.data,1,mean)
fNeglect.raw <- with(fug,data.frame(Q2.2_4,Q2.3_2,Q2.4_2,Q2.5_4,Q2.6_1,Q2.7_2,Q2.8_4,Q2.9_2,Q2.10_4,Q2.11_1,Q2.12_1,Q2.13_1)) # read Neglect items
fNeglect.data <- apply(fNeglect.raw,2,convert)
fNeglect <- apply(fNeglect.data,1,mean)
fVoice.raw <- with(fug,data.frame(Q2.2_2,Q2.3_10,Q2.4_3,Q2.5_1,Q2.6_4,Q2.7_1,Q2.8_3,Q2.9_3,Q2.10_2,Q2.11_2,Q2.12_4,Q2.13_3)) # read Voice items
fVoice.data <- apply(fVoice.raw,2,convert)
fVoice <- apply(fVoice.data,1,mean)
fLoyalty.raw <- with(fug,data.frame(Q2.2_3,Q2.3_1,Q2.4_4,Q2.5_3,Q2.6_3,Q2.7_3,Q2.8_1,Q2.9_1,Q2.10_3,Q2.11_3, Q2.12_2,Q2.13_2)) #read Loyalty items
fLoyalty.data <- apply(fLoyalty.raw,2,convert)
fLoyalty <- apply(fLoyalty.data,1,mean)

hop <- read.csv("./Data/Hoplock_data.csv", fileEncoding="latin1")
hop <- hop[-(1:2),]
hExclude <- convert(hop$Exclude)
hop <- hop[which(hExclude %in% include),]
hcond <- with(hop,convert(Q7.1))
hcond[hcond==1] <- 'H'
hcond[hcond==2] <- 'L'
hcond <- factor(hcond,levels=c('H','L'))
contrasts(hcond) <- contr.sum(2)
hExit.raw <- with(hop,data.frame(Q2.2_1,Q2.3_3,Q2.4_1,Q2.5_2,Q2.6_2,Q2.7_4,Q2.8_2,Q2.9_4,Q2.10_1,Q2.11_4,Q2.12_3,Q2.13_4)) # read Exit items
hExit.data <- apply(hExit.raw,2,convert)
hExit <- apply(hExit.data,1,mean)
hNeglect.raw <- with(hop,data.frame(Q2.2_4,Q2.3_2,Q2.4_2,Q2.5_4,Q2.6_1,Q2.7_2,Q2.8_4,Q2.9_2,Q2.10_4,Q2.11_1,Q2.12_1,Q2.13_1)) # read Neglect items
hNeglect.data <- apply(hNeglect.raw,2,convert)
hNeglect <- apply(hNeglect.data,1,mean)
hVoice.raw <- with(hop,data.frame(Q2.2_2,Q2.3_10,Q2.4_3,Q2.5_1,Q2.6_4,Q2.7_1,Q2.8_3,Q2.9_3,Q2.10_2,Q2.11_2,Q2.12_4,Q2.13_3)) # read Voice items
hVoice.data <- apply(hVoice.raw,2,convert)
hVoice <- apply(hVoice.data,1,mean)
hLoyalty.raw <- with(hop,data.frame(Q2.2_3,Q2.3_1,Q2.4_4,Q2.5_3,Q2.6_3,Q2.7_3,Q2.8_1,Q2.9_1,Q2.10_3,Q2.11_3, Q2.12_2,Q2.13_2)) #read Loyalty items
hLoyalty.data <- apply(hLoyalty.raw,2,convert)
hLoyalty <- apply(hLoyalty.data,1,mean)

sinclair <- read.csv("./Data/Sinclair_data.csv", fileEncoding="latin1")
sinclair <- sinclair[-(1:2),]
siExclude <- convert(sinclair$Exclude)
sinclair <- sinclair[which(siExclude %in% include),]
sicond <- with(sinclair,convert(Q7.1))
sicond[sicond==1] <- 'H'
sicond[sicond==2] <- 'L'
sicond <- factor(sicond,levels=c('H','L'))
contrasts(sicond) <- contr.sum(2)
siExit.raw <- with(sinclair,data.frame(Q2.2_1,Q2.3_3,Q2.4_1,Q2.5_2,Q2.6_2,Q2.7_4,Q2.8_2,Q2.9_4,Q2.10_1,Q2.11_4,Q2.12_3,Q2.13_4)) # read Exit items
siExit.data <- apply(siExit.raw,2,convert)
siExit <- apply(siExit.data,1,mean)
siNeglect.raw <- with(sinclair,data.frame(Q2.2_4,Q2.3_2,Q2.4_2,Q2.5_4,Q2.6_1,Q2.7_2,Q2.8_4,Q2.9_2,Q2.10_4,Q2.11_1,Q2.12_1,Q2.13_1)) # read Neglect items
siNeglect.data <- apply(siNeglect.raw,2,convert)
siNeglect <- apply(siNeglect.data,1,mean)
siVoice.raw <- with(sinclair,data.frame(Q2.2_2,Q2.3_10,Q2.4_3,Q2.5_1,Q2.6_4,Q2.7_1,Q2.8_3,Q2.9_3,Q2.10_2,Q2.11_2,Q2.12_4,Q2.13_3)) # read Voice items
siVoice.data <- apply(siVoice.raw,2,convert)
siVoice <- apply(siVoice.data,1,mean)
siLoyalty.raw <- with(sinclair,data.frame(Q2.2_3,Q2.3_1,Q2.4_4,Q2.5_3,Q2.6_3,Q2.7_3,Q2.8_1,Q2.9_1,Q2.10_3,Q2.11_3, Q2.12_2,Q2.13_2)) #read Loyalty items
siLoyalty.data <- apply(siLoyalty.raw,2,convert)
siLoyalty <- apply(siLoyalty.data,1,mean)

such <- read.csv("./Data/Sucharyna_data.csv", fileEncoding="latin1")
such <- such[-(1:2),]
suExclude <- convert(such$Exclude)
such <- such[which(suExclude %in% include),]
sucond <- with(such,convert(Q7.1))
sucond[sucond==1] <- 'H'
sucond[sucond==2] <- 'L'
sucond <- factor(sucond,levels=c('H','L'))
contrasts(sucond) <- contr.sum(2)
suExit.raw <- with(such,data.frame(Q2.2_1,Q2.3_3,Q2.4_1,Q2.5_2,Q2.6_2,Q2.7_4,Q2.8_2,Q2.9_4,Q2.10_1,Q2.11_4,Q2.12_3,Q2.13_4)) # read Exit items
suExit.data <- apply(suExit.raw,2,convert)
suExit <- apply(suExit.data,1,mean)
suNeglect.raw <- with(such,data.frame(Q2.2_4,Q2.3_2,Q2.4_2,Q2.5_4,Q2.6_1,Q2.7_2,Q2.8_4,Q2.9_2,Q2.10_4,Q2.11_1,Q2.12_1,Q2.13_1)) # read Neglect items
suNeglect.data <- apply(suNeglect.raw,2,convert)
suNeglect <- apply(suNeglect.data,1,mean)
suVoice.raw <- with(such,data.frame(Q2.2_2,Q2.3_10,Q2.4_3,Q2.5_1,Q2.6_4,Q2.7_1,Q2.8_3,Q2.9_3,Q2.10_2,Q2.11_2,Q2.12_4,Q2.13_3)) # read Voice items
suVoice.data <- apply(suVoice.raw,2,convert)
suVoice <- apply(suVoice.data,1,mean)
suLoyalty.raw <- with(such,data.frame(Q2.2_3,Q2.3_1,Q2.4_4,Q2.5_3,Q2.6_3,Q2.7_3,Q2.8_1,Q2.9_1,Q2.10_3,Q2.11_3, Q2.12_2,Q2.13_2)) #read Loyalty items
suLoyalty.data <- apply(suLoyalty.raw,2,convert)
suLoyalty <- apply(suLoyalty.data,1,mean)

tid <- read.csv("./Data/Tidwell_data.csv", fileEncoding="latin1")
tid <- tid[-(1:2),]
tExclude <- convert(tid$Exclude)
tid <- tid[which(tExclude %in% include),]
tcond <- with(tid,convert(Q7.1))
tcond[tcond==1] <- 'H'
tcond[tcond==2] <- 'L'
tcond <- factor(tcond,levels=c('H','L'))
contrasts(tcond) <- contr.sum(2)
tExit.raw <- with(tid,data.frame(Q2.2_1,Q2.3_3,Q2.4_1,Q2.5_2,Q2.6_2,Q2.7_4,Q2.8_2,Q2.9_4,Q2.10_1,Q2.11_4,Q2.12_3,Q2.13_4)) # read Exit items
tExit.data <- apply(tExit.raw,2,convert)
tExit <- apply(tExit.data,1,mean)
tNeglect.raw <- with(tid,data.frame(Q2.2_4,Q2.3_2,Q2.4_2,Q2.5_4,Q2.6_1,Q2.7_2,Q2.8_4,Q2.9_2,Q2.10_4,Q2.11_1,Q2.12_1,Q2.13_1)) # read Neglect items
tNeglect.data <- apply(tNeglect.raw,2,convert)
tNeglect <- apply(tNeglect.data,1,mean)
tVoice.raw <- with(tid,data.frame(Q2.2_2,Q2.3_10,Q2.4_3,Q2.5_1,Q2.6_4,Q2.7_1,Q2.8_3,Q2.9_3,Q2.10_2,Q2.11_2,Q2.12_4,Q2.13_3)) # read Voice items
tVoice.data <- apply(tVoice.raw,2,convert)
tVoice <- apply(tVoice.data,1,mean)
tLoyalty.raw <- with(tid,data.frame(Q2.2_3,Q2.3_1,Q2.4_4,Q2.5_3,Q2.6_3,Q2.7_3,Q2.8_1,Q2.9_1,Q2.10_3,Q2.11_3, Q2.12_2,Q2.13_2)) #read Loyalty items
tLoyalty.data <- apply(tLoyalty.raw,2,convert)
tLoyalty <- apply(tLoyalty.data,1,mean)

vranka <- read.csv("./Data/Vranka_data.csv", fileEncoding="latin1")
vranka <- vranka[-(1:2),]
vExclude <- convert(vranka$Exclude)
vranka <- vranka[which(vExclude %in% include),]
vcond <- with(vranka,convert(Q7.1))
vcond[vcond==1] <- 'H'
vcond[vcond==2] <- 'L'
vcond <- factor(vcond,levels=c('H','L'))
contrasts(vcond) <- contr.sum(2)
vExit.raw <- with(vranka,data.frame(Q2.2_1,Q2.3_3,Q2.4_1,Q2.5_2,Q2.6_2,Q2.7_4,Q2.8_2,Q2.9_4,Q2.10_1,Q2.11_4,Q2.12_3,Q2.13_4)) # read Exit items
vExit.data <- apply(vExit.raw,2,convert)
vExit <- apply(vExit.data,1,mean)
vNeglect.raw <- with(vranka,data.frame(Q2.2_4,Q2.3_2,Q2.4_2,Q2.5_4,Q2.6_1,Q2.7_2,Q2.8_4,Q2.9_2,Q2.10_4,Q2.11_1,Q2.12_1,Q2.13_1)) # read Neglect items
vNeglect.data <- apply(vNeglect.raw,2,convert)
vNeglect <- apply(vNeglect.data,1,mean)
vVoice.raw <- with(vranka,data.frame(Q2.2_2,Q2.3_10,Q2.4_3,Q2.5_1,Q2.6_4,Q2.7_1,Q2.8_3,Q2.9_3,Q2.10_2,Q2.11_2,Q2.12_4,Q2.13_3)) # read Voice items
vVoice.data <- apply(vVoice.raw,2,convert)
vVoice <- apply(vVoice.data,1,mean)
vLoyalty.raw <- with(vranka,data.frame(Q2.2_3,Q2.3_1,Q2.4_4,Q2.5_3,Q2.6_3,Q2.7_3,Q2.8_1,Q2.9_1,Q2.10_3,Q2.11_3, Q2.12_2,Q2.13_2)) #read Loyalty items
vLoyalty.data <- apply(vLoyalty.raw,2,convert)
vLoyalty <- apply(vLoyalty.data,1,mean)

yong <- read.csv("./Data/Yong_data.csv", fileEncoding="latin1")
yong <- yong[-(1:2),]
yExclude <- convert(yong$Exclude)
yong <- yong[which(yExclude %in% include),]
ycond <- with(yong,convert(Q7.1))
ycond[ycond==1] <- 'H'
ycond[ycond==2] <- 'L'
ycond <- factor(ycond,levels=c('H','L'))
contrasts(ycond) <- contr.sum(2)
yExit.raw <- with(yong,data.frame(Q2.2_1,Q2.3_3,Q2.4_1,Q2.5_2,Q2.6_2,Q2.7_4,Q2.8_2,Q2.9_4,Q2.10_1,Q2.11_4,Q2.12_3,Q2.13_4)) # read Exit items
yExit.data <- apply(yExit.raw,2,convert)
yExit <- apply(yExit.data,1,mean)
yNeglect.raw <- with(yong,data.frame(Q2.2_4,Q2.3_2,Q2.4_2,Q2.5_4,Q2.6_1,Q2.7_2,Q2.8_4,Q2.9_2,Q2.10_4,Q2.11_1,Q2.12_1,Q2.13_1)) # read Neglect items
yNeglect.data <- apply(yNeglect.raw,2,convert)
yNeglect <- apply(yNeglect.data,1,mean)
yVoice.raw <- with(yong,data.frame(Q2.2_2,Q2.3_10,Q2.4_3,Q2.5_1,Q2.6_4,Q2.7_1,Q2.8_3,Q2.9_3,Q2.10_2,Q2.11_2,Q2.12_4,Q2.13_3)) # read Voice items
yVoice.data <- apply(yVoice.raw,2,convert)
yVoice <- apply(yVoice.data,1,mean)
yLoyalty.raw <- with(yong,data.frame(Q2.2_3,Q2.3_1,Q2.4_4,Q2.5_3,Q2.6_3,Q2.7_3,Q2.8_1,Q2.9_1,Q2.10_3,Q2.11_3, Q2.12_2,Q2.13_2)) #read Loyalty items
yLoyalty.data <- apply(yLoyalty.raw,2,convert)
yLoyalty <- apply(yLoyalty.data,1,mean)

#Results
#Exit
exits = list(
  list(ayExit, aycond), list(bExit,  bcond), list(capExit, capcond), list(carExit, carcond), 
  list(carsExit, carscond), list(chExit, chcond), list(coExit, cocond), list(collinExit, collincond), 
  list(dExit,  dcond), list(fExit,  fcond), list(hExit,  hcond), list(siExit, sicond), 
  list(suExit, sucond), list(tExit,  tcond), list(vExit,  vcond), list(yExit,  ycond))

exit_es = do.call(rbind, lapply(exits, FUN=function(x){
    tmp = data.frame(y=x[[1]], cond=x[[2]]) %>% 
      filter(!is.na(y)) %>% group_by(cond) %>%
      summarize(mn = mean(y), vr = var(y), n=n()) %>% arrange(as.integer(cond == 'L'))
      d = (tmp$mn[1] - tmp$mn[2])/sqrt(sum((tmp$n-1) * tmp$vr)/(sum(tmp$n) - 2))
    return(data.frame(
      d = d, 
      vd = sum(tmp$n)/prod(tmp$n) + d^2/(2*sum(tmp$n)),
      ntrt = tmp$n[1],
      nctrl = tmp$n[2],
      t = d * sqrt(prod(tmp$n)/sum(tmp$n)),
      es = "d"
    ))
  })) %>%
  rbind(., data.frame(d = -.65/ sqrt(0.23^2 * (87 + 2)/4), 
                      vd = (-.65/ sqrt(0.23^2 * (87 + 2)/4))^2 / (2*89) + 4/89,
                      ntrt = 89/2,
                      nctrl = 89/2,
                      t = -.65/.23,
                      es = 'd'
                      ))
exit_es$site = c("Aykutoglu", "Bredow", "Caprariello", "Carcedo", 
                 "Carson", "Cheung", "Cobb", "Collins", "DiDonato", 
                 "Fuglestad", "Hoplock", "Sinclair", "Sucharyna", 
                 "Tidwell", "Vranka", "Yong", "original")
exit_es$experiment = "Exit"
exit_es$replicated = 0

# Original data
# diff = -.65, SE=.23, Df = 87, t = -.65/.23


#Neglect
# Original study
# Diff = -.42, SE = 0.199, Df = 87
orig_neglect = list(diff=-.42, se = .199, df = 87)

neglects = list(
  list(ayNeglect, aycond), list(bNeglect,  bcond), list(capNeglect, capcond), list(carNeglect, carcond), 
  list(carsNeglect, carscond), list(chNeglect, chcond), list(coNeglect, cocond), list(collinNeglect, collincond), 
  list(dNeglect,  dcond), list(fNeglect,  fcond), list(hNeglect,  hcond), list(siNeglect, sicond), 
  list(suNeglect, sucond), list(tNeglect,  tcond), list(vNeglect,  vcond), list(yNeglect,  ycond))

neglect_es = do.call(rbind, lapply(neglects, FUN=function(x){
  tmp = data.frame(y=x[[1]], cond=x[[2]]) %>% 
    filter(!is.na(y)) %>% group_by(cond) %>%
    summarize(mn = mean(y), vr = var(y), n=n()) %>% arrange(as.integer(cond == 'L'))
  d = (tmp$mn[1] - tmp$mn[2])/sqrt(sum((tmp$n-1) * tmp$vr)/(sum(tmp$n) - 2))
  return(data.frame(
    d = d, 
    vd = sum(tmp$n)/prod(tmp$n) + d^2/(2*sum(tmp$n)),
    ntrt = tmp$n[1],
    nctrl = tmp$n[2],
    t = d * sqrt(prod(tmp$n)/sum(tmp$n)),
    es = "d"
  ))
})) %>%
  rbind(., data.frame(d = orig_neglect$diff/ sqrt(orig_neglect$se^2 * (orig_neglect$df + 2)/4), 
                      vd = (orig_neglect$diff/ sqrt(orig_neglect$se^2 * (orig_neglect$df + 2)/4))^2 / (2 * orig_neglect$df + 4) + 4/(orig_neglect$df + 2),
                      ntrt = (orig_neglect$df + 2)/2,
                      nctrl = (orig_neglect$df + 2)/2,
                      t = orig_neglect$diff/orig_neglect$se,
                      es = 'd'
  ))
neglect_es$site = c("Aykutoglu", "Bredow", "Caprariello", "Carcedo", 
                    "Carson", "Cheung", "Cobb", "Collins", "DiDonato", 
                    "Fuglestad", "Hoplock", "Sinclair", "Sucharyna", 
                    "Tidwell", "Vranka", "Yong", "original")
neglect_es$experiment = "Neglect"
neglect_es$replicated = 0

forest(neglect_es$d, vi=neglect_es$vd)


#Voice
#original study
# Diff = .46, SE = 0.2857, Df = 87
orig_voice = list(diff=.46, se = .2857, df = 87)
voices = list(
  list(ayVoice, aycond), list(bVoice,  bcond), list(capVoice, capcond), list(carVoice, carcond), 
  list(carsVoice, carscond), list(chVoice, chcond), list(coVoice, cocond), list(collinVoice, collincond), 
  list(dVoice,  dcond), list(fVoice,  fcond), list(hVoice,  hcond), list(siVoice, sicond), 
  list(suVoice, sucond), list(tVoice,  tcond), list(vVoice,  vcond), list(yVoice,  ycond))

voice_es = do.call(rbind, lapply(voices, FUN=function(x){
  tmp = data.frame(y=x[[1]], cond=x[[2]]) %>% 
    filter(!is.na(y)) %>% group_by(cond) %>%
    summarize(mn = mean(y), vr = var(y), n=n()) %>% arrange(as.integer(cond == 'L'))
  d = (tmp$mn[1] - tmp$mn[2])/sqrt(sum((tmp$n-1) * tmp$vr)/(sum(tmp$n) - 2))
  return(data.frame(
    d = d, 
    vd = sum(tmp$n)/prod(tmp$n) + d^2/(2*sum(tmp$n)),
    ntrt = tmp$n[1],
    nctrl = tmp$n[2],
    t = d * sqrt(prod(tmp$n)/sum(tmp$n)),
    es = "d"
  ))
})) %>%
  rbind(., data.frame(d = orig_voice$diff/ sqrt(orig_voice$se^2 * (orig_voice$df + 2)/4), 
                      vd = (orig_voice$diff/ sqrt(orig_voice$se^2 * (orig_voice$df + 2)/4))^2 / (2 * orig_voice$df + 4) + 4/(orig_voice$df + 2),
                      ntrt = (orig_voice$df + 2)/2,
                      nctrl = (orig_voice$df + 2)/2,
                      t = orig_voice$diff/orig_voice$se,
                      es = 'd'
  ))
voice_es$site = c("Aykutoglu", "Bredow", "Caprariello", "Carcedo", 
                  "Carson", "Cheung", "Cobb", "Collins", "DiDonato", 
                  "Fuglestad", "Hoplock", "Sinclair", "Sucharyna", 
                  "Tidwell", "Vranka", "Yong", "original")
voice_es$experiment = "Voice"
voice_es$replicated = 1

forest(voice_es$d, vi=voice_es$vd)


#Loyalty
#original study
# Diff = .29, SE = 0.26, Df = 87
orig_loyalty = list(diff=.26, se = .26, df = 87)
loyals = list(
  list(ayLoyalty, aycond), list(bLoyalty,  bcond), list(capLoyalty, capcond), list(carLoyalty, carcond), 
  list(carsLoyalty, carscond), list(chLoyalty, chcond), list(coLoyalty, cocond), list(collinLoyalty, collincond), 
  list(dLoyalty,  dcond), list(fLoyalty,  fcond), list(hLoyalty,  hcond), list(siLoyalty, sicond), 
  list(suLoyalty, sucond), list(tLoyalty,  tcond), list(vLoyalty,  vcond), list(yLoyalty,  ycond))

loyalty_es = do.call(rbind, lapply(loyals, FUN=function(x){
  tmp = data.frame(y=x[[1]], cond=x[[2]]) %>% 
    filter(!is.na(y)) %>% group_by(cond) %>%
    summarize(mn = mean(y), vr = var(y), n=n()) %>% arrange(as.integer(cond == 'L'))
  d = (tmp$mn[1] - tmp$mn[2])/sqrt(sum((tmp$n-1) * tmp$vr)/(sum(tmp$n) - 2))
  return(data.frame(
    d = d, 
    vd = sum(tmp$n)/prod(tmp$n) + d^2/(2*sum(tmp$n)),
    ntrt = tmp$n[1],
    nctrl = tmp$n[2],
    t = d * sqrt(prod(tmp$n)/sum(tmp$n)),
    es = "d"
  ))
})) %>%
  rbind(., data.frame(d = orig_loyalty$diff/ sqrt(orig_loyalty$se^2 * (orig_loyalty$df + 2)/4), 
                      vd = (orig_loyalty$diff/ sqrt(orig_loyalty$se^2 * (orig_loyalty$df + 2)/4))^2 / (2 * orig_loyalty$df + 4) + 4/(orig_loyalty$df + 2),
                      ntrt = (orig_loyalty$df + 2)/2,
                      nctrl = (orig_loyalty$df + 2)/2,
                      t = orig_loyalty$diff/orig_loyalty$se,
                      es = 'd'
  ))
loyalty_es$site = c("Aykutoglu", "Bredow", "Caprariello", "Carcedo", 
                    "Carson", "Cheung", "Cobb", "Collins", "DiDonato", 
                    "Fuglestad", "Hoplock", "Sinclair", "Sucharyna", 
                    "Tidwell", "Vranka", "Yong", "original")
loyalty_es$experiment = "Loyalty"
loyalty_es$replicated = 1

forest(loyalty_es$d, vi=loyalty_es$vd)


# combine dfs and write
df <- bind_rows(exit_es, neglect_es, voice_es, loyalty_es)

write.csv(df,  "../../rrr_cheung.csv", row.names=F)


