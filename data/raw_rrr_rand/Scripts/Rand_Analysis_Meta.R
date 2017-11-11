source('Rand_Analysis.R') # load analysis script from source location
if(suppressWarnings(!require(metafor,quietly=TRUE))) install.packages('metafor') # install metafor package if necessary
library(metafor) # load metafor package

#======#
# Data #
#======#
Original <- read.table('Original_data.txt',header=TRUE)[,-c(1,5,8,9)]

csv <- list.files(pattern="*.csv") # names of all csv files in folder
authors <- sapply(strsplit(csv,split='_',fixed=TRUE),function(x){x[1]}) # authors' names
k <- length(authors) # number of labs
data.all <- list()
for(i in 1:k){
  data <- read.csv(csv[i],header=FALSE) # read csv file
  if(as.character(data[1,1])=='startDate'){ # if data in new format
    data.all[[i]] <- data[,c(1:10,14:22,11,23:27,12,28:30,13,31:33,13,34:87)] # make compatible
  } else data.all[[i]] <- data # otherwise, no change
}
names(data.all) <- authors # assign author name to each dataset in list
Exclusions <- c('None','Experienced','Non-Compliant','Non-Comprehending','Experienced, Non-Compliant, or Non-Comprehending')

tab <- as.data.frame(t(sapply(1:k,function(i){Rand.DS(data.all[[i]])})),row.names=authors) # table of descriptive stats
colnames(tab) <- c('N_Total_TP','N_Women_TP','M_Age_TP','SD_Age_TP','N_Personal_TP','N_Group_TP','M_Trust1_TP','SD_Trust1_TP','M_Trust2_TP','SD_Trust2_TP',
                   'N_Trust3_TP','N_Trust4_TP','N_Trust5_TP','N_Trust6_TP','N_Exp1_TP','M_Exp2_TP','SD_Exp2_TP','M_Exp3_TP','SD_Exp3_TP','M_Exp4_TP','SD_Exp4_TP',
                   'M_Exp5_TP','SD_Exp5_TP','M_KOP_TP','SD_KOP_TP','M_HI_TP','SD_HI_TP','M_VI_TP','SD_VI_TP','M_HC_TP','SD_HC_TP','M_VC_TP','SD_VC_TP',
                   'N_Total_FD','N_Women_FD','M_Age_FD','SD_Age_FD','N_Personal_FD','N_Group_FD','M_Trust1_FD','SD_Trust1_FD','M_Trust2_FD','SD_Trust2_FD',
                   'N_Trust3_FD','N_Trust4_FD','N_Trust5_FD','N_Trust6_FD','N_Exp1_FD','M_Exp2_FD','SD_Exp2_FD','M_Exp3_FD','SD_Exp3_FD','M_Exp4_FD','SD_Exp4_FD',
                   'M_Exp5_FD','SD_Exp5_FD','M_KOP_FD','SD_KOP_FD','M_HI_FD','SD_HI_FD','M_VI_FD','SD_VI_FD','M_HC_FD','SD_HC_FD','M_VC_FD','SD_VC_FD',
                   'N_All_TP','M_Time_All_TP','SD_Time_All_TP','M_Cont_All_TP','SD_Cont_All_TP','N_All_FD','M_Time_All_FD','SD_Time_All_FD','M_Cont_All_FD','SD_Cont_All_FD',
                   'N_E1_TP','M_Time_E1_TP','SD_Time_E1_TP','M_Cont_E1_TP','SD_Cont_E1_TP','N_E1_FD','M_Time_E1_FD','SD_Time_E1_FD','M_Cont_E1_FD','SD_Cont_E1_FD',
                   'N_E2_TP','M_Time_E2_TP','SD_Time_E2_TP','M_Cont_E2_TP','SD_Cont_E2_TP','N_E2_FD','M_Time_E2_FD','SD_Time_E2_FD','M_Cont_E2_FD','SD_Cont_E2_FD',
                   'N_E3_TP','M_Time_E3_TP','SD_Time_E3_TP','M_Cont_E3_TP','SD_Cont_E3_TP','N_E3_FD','M_Time_E3_FD','SD_Time_E3_FD','M_Cont_E3_FD','SD_Cont_E3_FD',
                   'N_E123_TP','M_Time_E123_TP','SD_Time_E123_TP','M_Cont_E123_TP','SD_Cont_E123_TP','N_E123_FD','M_Time_E123_FD','SD_Time_E123_FD','M_Cont_E123_FD','SD_Cont_E123_FD')
write.csv(tab,'Table.csv')

#======#
# Main #
#======#

#----------#
# Original #
#----------#
sink('Original_Main.txt')
for(i in c(1,3,4)){
  cat('\nExclude: ',Exclusions[i],'\n',rep('=',nchar(Exclusions[i])+9),'\n\n',sep='')
  if(i==1) O.data <- Original
  if(i==3) O.data <- Original[Original$Disobeyed==0,]
  if(i==4) O.data <- Original[Original$Failed_Comprehension==0,]
  O.Contribution <- O.data$Contribution/4
  O.Condition <- factor(ifelse(O.data$Condition=='Time_delay','FD','TP'),levels=c('TP','FD'))
  contrasts(O.Condition) <- contr.sum(2)
  ttest(O.Contribution,O.Condition)
}
sink()

#------#
# Labs #
#------#
for(j in 1:k){
  out <- paste(authors[j],'_Main.txt',sep='') # start writing named output file
  sink(out)
  for(i in 1:5){
    cat('\nExclude: ',Exclusions[i],'\n',rep('=',nchar(Exclusions[i])+9),'\n\n',sep='')
    Rand.Analysis(data.all[[j]],type='main',exclude=i-1)
  }
  sink()
}

#------#
# Meta #
#------#

#RMA w/ Original
authors.0 <- c('Original',authors)
Main.0 <- HL <- vector('list',5)
temp <- data.frame(matrix(0,k+1,4))
dimnames(temp) <- list(authors.0,c('D','SE','TP','FD'))
for(i in c(1,3,4)){
  if(i==1) O.data <- Original
  if(i==3) O.data <- Original[Original$Disobeyed==0,]
  if(i==4) O.data <- Original[Original$Failed_Comprehension==0,]
  for(j in 1:(k+1)){
    if(j==1) {
      O.Contribution <- O.data$Contribution/4
      O.Condition <- factor(ifelse(O.data$Condition=='Time_delay','FD','TP'),levels=c('TP','FD'))
      contrasts(O.Condition) <- contr.sum(2)
      temp[1,] <- ttest(O.Contribution,O.Condition)
    } else temp[j,] <- Rand.Analysis(data.all[[j-1]],type='main',exclude=i-1)[,-5]
  }
  jn <- which(!is.na(temp$D))
  temp <- temp[jn,]
  Main.0[[i]] <- rma(yi=temp$D,sei=temp$SE,slab=authors.0[jn],knha=TRUE,control=list(stepadj=.5))
  HL[[i]] <- temp[,3:4]
}

# RMA w/o Original #
Main <- Main.T <- dat <- xa0 <- vector('list',5)
temp <- data.frame(matrix(0,k,5))
dimnames(temp) <- list(authors,c('D','SE','TP','FD','Trust'))
for(i in 1:5){
  for(j in 1:k){
    temp[j,] <- Rand.Analysis(data.all[[j]],type='main',exclude=i-1)
  }
  jn <- which(!is.na(temp$D))
  temp <- temp[jn,]
  Main[[i]] <- rma(yi=temp$D,sei=temp$SE,slab=authors[jn],knha=TRUE,control=list(stepadj=.5))
  Main.T[[i]] <- rma(yi=temp$D,sei=temp$SE,mods=temp$Trust,slab=authors[jn],knha=TRUE,control=list(stepadj=.5))
  if(i %in% c(2,5)) HL[[i]] <- temp[,3:4]
  dat[[i]] <- temp
  if(length(jn)<k){
    xa0[[i]] <- authors[-jn]
  } else xa0[[i]] <- NA
}
sink('RMA_Main.txt')
cat('META-ANALYSES OF EFFECT SIZE \n\n')
cat('Effect Size: Difference in Contribution Means between Conditions\n\n')
for(i in 1:5){
  cat('__________________________________________________________________________________\n')
  cat('\nExclude: ',Exclusions[i],'\n',rep('=',nchar(Exclusions[i])+9),'\n',sep='')
  if(!is.na(xa0[[i]])[1]) cat('*Studies Omitted For Insufficient Data:',paste(xa0[[i]],collapse=', '),'\n')
  print(Main[[i]],signif.legend=FALSE,signif.stars=FALSE)
}
sink()

# Forest Plot #
par(mar=c(6.1,5.1,4.1,3.1))
forest.main <- vector('list',5)

H <- round(HL[[1]]$TP,2)
L <- round(HL[[1]]$FD,2)
defaults <- forest(Main.0[[1]],cex=.7)
kk <- max(defaults$rows)-1
xp1 <- defaults$xlim[1]+.175*(defaults$xlim[2]-defaults$xlim[1])
xp2 <- defaults$xlim[1]+.255*(defaults$xlim[2]-defaults$xlim[1])
DV.forest <- forest(Main.0[[1]],ilab=cbind(H[!is.na(H)],L[!is.na(L)]),ilab.xpos=c(xp1,xp2),xlab='Effect Size (Difference in Means)',addfit=FALSE,
                    ylim=c(-10,kk+3),rows=c(kk-1,seq(kk-3,-2,-1)),cex=.7)
text(defaults$xlim[1],kk,adj=c(0,0),'STUDY',font=2,cex=.7) # column name for labs
text(defaults$xlim[1]+.155*(defaults$xlim[2]-defaults$xlim[1]),kk+.35,adj=c(0,0),'MEAN CONTRIBUTION',font=2,cex=.7) # column header
text(defaults$xlim[1]+.1525*(defaults$xlim[2]-defaults$xlim[1]),kk-.35,adj=c(0,0),'Pressure',font=2,cex=.7) # column name for TP means
text(defaults$xlim[1]+.2425*(defaults$xlim[2]-defaults$xlim[1]),kk-.35,adj=c(0,0),'Delay',font=2,cex=.7) # column name for FD means
text(defaults$xlim[1]+.9225*(defaults$xlim[2]-defaults$xlim[1]),kk+.35,adj=c(0,0),'DIFFERENCE',font=2,cex=.7) # column name for confidence intervals
text(defaults$xlim[1]+.8575*(defaults$xlim[2]-defaults$xlim[1]),kk-.35,adj=c(0,0),'Pressure - Delay [95% CI]',font=2,cex=.7) # column name for confidence intervals
if(!is.na(xa0[[1]])[1]) text(defaults$xlim[1],kk+1.25,adj=c(0,0),paste('*Studies omitted for insufficient data:',paste(xa0[[i]],collapse=', ')),cex=.7)
abline(h=kk-2,lty=2)
abline(h=-3)
text(defaults$xlim[1],-4,adj=c(0,0.3),'SUMMARY (Random Effects)',font=2,cex=.7)
addpoly(Main[[1]],row=-5,mlab=paste('Exclude:',Exclusions[1]),font=2,cex=.7)
for(i in 2:5) addpoly(Main[[i]],row=-(i+4.5),mlab=paste('Exclude:',Exclusions[i]),cex=.7)
forest.main[[1]] <- recordPlot()

for(i in 2:5){
  H <- round(HL[[i]]$TP,2)
  L <- round(HL[[i]]$FD,2)
  if(i %in% c(3,4)){
    defaults <- forest(Main.0[[i]],cex=.8)
    kk <- max(defaults$rows)-1
    xp1 <- defaults$xlim[1]+.175*(defaults$xlim[2]-defaults$xlim[1]) # x-coordinate for high means column
    xp2 <- defaults$xlim[1]+.255*(defaults$xlim[2]-defaults$xlim[1])
    DV.forest <- forest(Main.0[[i]],ilab=cbind(H[!is.na(H)],L[!is.na(L)]),ilab.xpos=c(xp1,xp2),xlab='Effect Size (Difference in Means)',addfit=FALSE,
                        ylim=c(-5,kk+3),rows=c(kk-1,seq(kk-3,-2,-1)),cex=.8)
    text(defaults$xlim[1],kk,adj=c(0,0),'STUDY',font=2,cex=.8)
    text(defaults$xlim[1]+.1355*(defaults$xlim[2]-defaults$xlim[1]),kk+.35,adj=c(0,0),'MEAN CONTRIBUTION',font=2,cex=.8) # column header
    text(defaults$xlim[1]+.1455*(defaults$xlim[2]-defaults$xlim[1]),kk-.35,adj=c(0,0),'Pressure',font=2,cex=.8) # column name for TP means
    text(defaults$xlim[1]+.2375*(defaults$xlim[2]-defaults$xlim[1]),kk-.35,adj=c(0,0),'Delay',font=2,cex=.8) # column name for FD means
    text(defaults$xlim[1]+.9050*(defaults$xlim[2]-defaults$xlim[1]),kk+.35,adj=c(0,0),'DIFFERENCE',font=2,cex=.8) # column name for confidence intervals
    text(defaults$xlim[1]+.8250*(defaults$xlim[2]-defaults$xlim[1]),kk-.35,adj=c(0,0),'Pressure - Delay [95% CI]',font=2,cex=.8) # column name for confidence intervals
    text(mean(defaults$xlim),kk+3,paste('Exclude:',Exclusions[i]),font=2,cex=1)
    if(!is.na(xa0[[i]])[1]) text(defaults$xlim[1],kk+1.25,adj=c(0,0),paste('*Studies omitted for insufficient data:',paste(xa0[[i]],collapse=', ')),cex=.7)
    abline(h=kk-2,lty=2)
    abline(h=-3)
    addpoly(Main[[i]],row=-4,mlab='SUMMARY (Random Effects)',font=2,cex=.8)
  } else {
    defaults <- forest(Main[[i]],cex=.8)
    kk <- max(defaults$rows)
    xp1 <- defaults$xlim[1]+.175*(defaults$xlim[2]-defaults$xlim[1]) # x-coordinate for high means column
    xp2 <- defaults$xlim[1]+.255*(defaults$xlim[2]-defaults$xlim[1])
    DV.forest <- forest(Main[[i]],ilab=cbind(H[!is.na(H)],L[!is.na(L)]),ilab.xpos=c(xp1,xp2),xlab='Effect Size (Difference in Means)',addfit=FALSE,
                        ylim=c(-3,kk+3),rows=seq(kk-1,0,-1),cex=.8)
    text(defaults$xlim[1],kk,adj=c(0,0),'STUDY',font=2,cex=.8) # column name for labs
    text(defaults$xlim[1]+.1355*(defaults$xlim[2]-defaults$xlim[1]),kk+.35,adj=c(0,0),'MEAN CONTRIBUTION',font=2,cex=.8) # column header
    text(defaults$xlim[1]+.1455*(defaults$xlim[2]-defaults$xlim[1]),kk-.35,adj=c(0,0),'Pressure',font=2,cex=.8) # column name for TP means
    text(defaults$xlim[1]+.2375*(defaults$xlim[2]-defaults$xlim[1]),kk-.35,adj=c(0,0),'Delay',font=2,cex=.8) # column name for FD means
    text(defaults$xlim[1]+.9050*(defaults$xlim[2]-defaults$xlim[1]),kk+.35,adj=c(0,0),'DIFFERENCE',font=2,cex=.8) # column name for confidence intervals
    text(defaults$xlim[1]+.8250*(defaults$xlim[2]-defaults$xlim[1]),kk-.35,adj=c(0,0),'Pressure - Delay [95% CI]',font=2,cex=.8) # column name for confidence intervals
    text(mean(defaults$xlim),kk+3,paste('Exclude:',Exclusions[i]),font=2,cex=1)
    if(!is.na(xa0[[i]])[1]) text(defaults$xlim[1],kk+1.25,adj=c(0,0),paste('*Studies omitted for insufficient data:',paste(xa0[[i]],collapse=', ')),cex=.7)
    abline(h=-1)
    addpoly(Main[[i]],row=-2,mlab='SUMMARY (Random Effects)',font=2,cex=.8)
  }
  forest.main[[i]] <- recordPlot()
}
pdf('Plots_Main.pdf',11,8.5) # start saving plots to pdf file
sapply(forest.main,replayPlot) # generate forest plot for each analysis
graphics.off() # end saving plots

#============#
# Moderators #
#============#

modnames <- rep(c('Trust','Gender','Age','Subject Pool Experience','Paid Experience','MTurk Experience','Deception Experience',
                  'Knowledge of Other Participants','Horizontal-Individualism','Vertical-Individualism','Horizontal-Collectivism',
                  'Vertical-Collectivism'),rep(2,12))

#----------#
# Original #
#----------#
sink('Original_Moderator_Exclude_None.txt')
O.Contribution <- Original$Contribution/4
O.Condition <- factor(ifelse(Original$Condition=='Time_delay','FD','TP'),levels=c('TP','FD'))
contrasts(O.Condition) <- contr.sum(2)
Gender <- factor(Original$Female+1)
contrasts(Gender) <- contr.sum(2)
Gender.out <- ancova(O.Contribution,O.Condition,Gender)$Interaction
Age <- Original$Age
Age.out <- ancova(O.Contribution,O.Condition,Age)$Interaction
sink()

#------#
# Labs #
#------#
for(i in c(1,5)){
  if(i==1) Exc <- 'None'
  if(i==5) Exc <- 'All'
  for(j in 1:k){
    out <- paste(authors[j],'_Moderator_Exclude_',Exc,'.txt',sep='') # start writing named output file
    sink(out)
    Rand.Analysis(data.all[[j]],type='moderator',exclude=i-1)
    sink()
  }
}

#------#
# Meta #
#------#
Mod0 <- Moderator0 <- xa0 <- vector('list',24)
for(i in 1:24){
  if(i %in% c(3,5)){
    Mod0[[i]] <- data.frame(matrix(0,k+1,4))
  } else{
    Mod0[[i]] <- data.frame(matrix(0,k,4))
  }
  names(Mod0[[i]]) <- c('D','SE','TP','FD')
}
Mod0[[3]][1,] <- Gender.out
Mod0[[5]][1,] <- Age.out
for(j in 1:k){
  out1 <- Rand.Analysis(data.all[[j]],type='moderator',exclude=0)
  out2 <- Rand.Analysis(data.all[[j]],type='moderator',exclude=4)
  for(l in 1:12){
    if(l %in% 2:3){
      Mod0[[l*2-1]][j+1,] <- out1[[l]]$Interaction
    } else{
      Mod0[[l*2-1]][j,] <- out1[[l]]$Interaction
    }
    Mod0[[l*2]][j,] <- out2[[l]]$Interaction
    names(Mod0)[l*2-1] <- paste(names(out1)[l],'.Exclude.None',sep='')
    names(Mod0)[l*2] <- paste(names(out2)[l],'.Exclude.All',sep='')
  }
}
for(i in 1:24){
  if(i %in% c(3,5)){
    temp <- Mod0[[i]][-1,]
  } else{
    temp <- Mod0[[i]]
  }
  j <- which(!is.na(temp$D))
  temp <- temp[j,]
  Moderator0[[i]] <- rma(yi=temp$D,sei=temp$SE,slab=authors[j],knha=TRUE,control=list(stepadj=.5))
  if(length(j)<k){
    xa0[[i]] <- authors[-j]
  } else xa0[[i]] <- NA
}
names(Moderator0) <- names(xa0) <- names(Mod0)

Moderator00 <- vector('list',2)
for(i in c(3,5)){
  temp <- Mod0[[i]]
  j <- which(!is.na(temp$D))
  temp <- temp[j,]
  Moderator00[[(i-1)/2]] <- rma(yi=temp$D,sei=temp$SE,slab=authors.0,knha=TRUE,control=list(stepadj=.5))
}

# RMA #
sink('RMA_Moderator.txt')
cat('META-REGRESSION OF EFFECT SIZE\n\n')
cat('Effect Size: Difference in Contribution Means between Conditions\n')
cat('Moderator: Average Trust of Study\n\n')
for(i in c(1,5)){
  cat('__________________________________________________________________________________\n')
  if(i==1){
    cat('\nNo Exclusions\n')
    cat('=============\n')
  } else{
    cat('\nAll Exclusions (Experienced, Non-Compliant, or Non-Comprehending)\n')
    cat('=================================================================\n')
    cat('*Studies Omitted For Insufficient Data:',paste(authors[-jn],collapse=', '),'\n')
  }
  print(Main.T[[i]],signif.legend=FALSE,signif.stars=FALSE)
}
cat('__________________________________________________________________________________\n')
cat('__________________________________________________________________________________\n\n\n')
cat('META-ANALYSES OF MODERATION EFFECT\n\n')
cat('For Continuous Moderators: Difference in Moderator Slopes between Conditions \n')
cat('                           (i.e., degree of moderator influence on effect size)\n')
cat('For Dichotomous Moderators: Difference in Effect Sizes between Moderator Levels \n\n')
for(i in 1:24){
  if(i %in% seq(1,23,2)){
    cat('__________________________________________________________________________________\n')
    cat('\n[',modnames[i],']\n\n',sep='')
    cat('No Exclusions\n')
    cat('=============\n')
  } else{
    cat('All Exclusions (Experienced, Non-Compliant, or Non-Comprehending)\n')
    cat('=================================================================\n')
  }
  if(!is.na(xa0[[i]])[1]) cat('*Studies Omitted For Insufficient Data:',paste(xa0[[i]],collapse=', '),'\n')
  print(Moderator0[[i]],signif.legend=FALSE,signif.stars=FALSE)
  if(i %in% seq(1,23,2)) cat('.............................................................................\n\n')
}
sink()

# Plots #
par(mar=c(6.1,5.1,6.1,3.1))
plots.mod <- vector('list',26)
ii <- 0

for(i in c(1,5)){
  ii <- ii+1
  dT <- dat[[i]]$Trust
  dD <- dat[[i]]$D
  sdT <- seq(min(dT)-.5,max(dT)+.5,.1)
  preds <- predict(Main.T[[i]],newmods=sdT)
  wi <- 1/dat[[i]]$SE
  size <- 3+4*(wi-min(wi))/(max(wi)-min(wi))
  yl <- min(c(preds$ci.lb,dD))
  yu <- max(c(preds$ci.ub,dD))
  if(i==1){
    subex <- 'No Exclusions'
    auths <- authors
  } 
  if(i==5){
    subex <- 'All Exclusions (Experienced, Non-Compliant, or Non-Comprehending)'
    auths <- authors[jn]
  } 
  plot(dT,dD,pch=21,col=NULL,bg='grey',cex=size,xlab="Average Trust",ylab="Effect Size",ylim=c(yl,yu),
       main=paste('Meta-Regression'))
  mtext(subex,font=2,padj=-2.2)
  if(length(auths)<length(authors)){
    mtext(paste('*Studies Omitted For Insufficient Data:',paste(authors[-jn],collapse=', ')),cex=.7,adj=0,padj=-.5)
  } 
  lines(sdT,preds$pred)
  lines(sdT,preds$ci.lb,lty="dashed")
  lines(sdT,preds$ci.ub,lty="dashed")
  text(dT,dD,auths,cex=.75)
  legend('top',c(paste('slope = ',round(Main.T[[i]]$b[2],2)),paste('p-value = ',round(Main.T[[i]]$pval[2],2))),bty='n',cex=.9)
  plots.mod[[ii]] <- recordPlot()
}

par(mar=c(6.1,5.1,4.1,3.1))

cont <- c(1,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1)
xl <- ifelse(cont==1,'Moderation Effect (Difference in Slopes)','Moderation Effect (Difference in Effect Sizes)')
cn <- ifelse(cont==1,'MODERATOR SLOPE','        EFFECT SIZE')
cn1 <- ifelse(cont==1,'Pressure','    Yes')
cn1[c(3,4)] <- '   Male'
cn2 <- ifelse(cont==1,'  Delay','    No')
cn2[c(3,4)] <- 'Female'
cn3 <- ifelse(cont==1,'Pressure - Delay [95% CI]','              Yes - No [95% CI]')
cn3[c(3,4)] <- '     Male - Female [95% CI]'
excl <- rep(c('No Exclusions','All Exclusions (Experienced, Non-Compliant, or Non-Comprehending)'),12)
for(i in 1:24){
  H <- round(Mod0[[i]]$TP,2)
  L <- round(Mod0[[i]]$FD,2)
  if(i %in% c(3,5)){
    defaults <- forest(Moderator00[[(i-1)/2]],cex=.8)
    kk <- max(defaults$rows)-1
    xp1 <- defaults$xlim[1]+.175*(defaults$xlim[2]-defaults$xlim[1])
    xp2 <- defaults$xlim[1]+.255*(defaults$xlim[2]-defaults$xlim[1])
    DV.forest <- forest(Moderator00[[(i-1)/2]],ilab=cbind(H[!is.na(H)],L[!is.na(L)]),ilab.xpos=c(xp1,xp2),xlab=xl[i],addfit=FALSE,
                        ylim=c(-5,kk+3),rows=c(kk-1,seq(kk-3,-2,-1)),cex=.8)
    text(defaults$xlim[1],kk,adj=c(0,0),'STUDY',font=2,cex=.8)
    text(defaults$xlim[1]+.1425*(defaults$xlim[2]-defaults$xlim[1]),kk+.35,adj=c(0,0),cn[i],font=2,cex=.8) # column header
    text(defaults$xlim[1]+.1475*(defaults$xlim[2]-defaults$xlim[1]),kk-.35,adj=c(0,0),cn1[i],font=2,cex=.8) # column name for TP means
    text(defaults$xlim[1]+.2300*(defaults$xlim[2]-defaults$xlim[1]),kk-.35,adj=c(0,0),cn2[i],font=2,cex=.8) # column name for FD means
    text(defaults$xlim[1]+.9075*(defaults$xlim[2]-defaults$xlim[1]),kk+.35,adj=c(0,0),'DIFFERENCE',font=2,cex=.8) # column name for confidence intervals
    text(defaults$xlim[1]+.8250*(defaults$xlim[2]-defaults$xlim[1]),kk-.35,adj=c(0,0),cn3[i],font=2,cex=.8) # column name for confidence intervals
    text(mean(defaults$xlim),kk+3.5,paste(modnames[i],'Moderator'),font=2,cex=1.25) # title indicating variable of interest
    text(mean(defaults$xlim),kk+2.5,excl[i],font=2,cex=1)
    if(!is.na(xa0[[i]])[1]) text(defaults$xlim[1],kk+1.25,adj=c(0,0),paste('*Studies omitted for insufficient data:',paste(xa0[[i]],collapse=', ')),cex=.7)
    abline(h=-3,lty=2)
    abline(h=kk-2,lty=2)
    addpoly(Moderator0[[i]],row=-4,mlab='SUMMARY (Random Effects)',font=2,cex=.8)
  } else{
    defaults <- forest(Moderator0[[i]],cex=.8)
    kk <- max(defaults$rows)
    xp1 <- defaults$xlim[1]+.175*(defaults$xlim[2]-defaults$xlim[1])
    xp2 <- defaults$xlim[1]+.255*(defaults$xlim[2]-defaults$xlim[1])
    DV.forest <- forest(Moderator0[[i]],ilab=cbind(H[!is.na(H)],L[!is.na(L)]),ilab.xpos=c(xp1,xp2),xlab=xl[i],addfit=FALSE,
                        ylim=c(-3,kk+3),rows=seq(kk-1,0,-1),cex=.8)
    text(defaults$xlim[1],kk,adj=c(0,0),'STUDY',font=2,cex=.8)
    text(defaults$xlim[1]+.1425*(defaults$xlim[2]-defaults$xlim[1]),kk+.35,adj=c(0,0),cn[i],font=2,cex=.8) # column header
    text(defaults$xlim[1]+.1475*(defaults$xlim[2]-defaults$xlim[1]),kk-.35,adj=c(0,0),cn1[i],font=2,cex=.8) # column name for TP means
    text(defaults$xlim[1]+.2300*(defaults$xlim[2]-defaults$xlim[1]),kk-.35,adj=c(0,0),cn2[i],font=2,cex=.8) # column name for FD means
    text(defaults$xlim[1]+.9075*(defaults$xlim[2]-defaults$xlim[1]),kk+.35,adj=c(0,0),'DIFFERENCE',font=2,cex=.8) # column name for confidence intervals
    text(defaults$xlim[1]+.8250*(defaults$xlim[2]-defaults$xlim[1]),kk-.35,adj=c(0,0),cn3[i],font=2,cex=.8) # column name for confidence intervals
    if(kk<10) kka <- .25 else kka <- .5
    text(mean(defaults$xlim),kk+2.85+kka,paste(modnames[i],'Moderator'),font=2,cex=1.25) # title indicating variable of interest
    text(mean(defaults$xlim),kk+1.85+kka,excl[i],font=2,cex=1)
    if(!is.na(xa0[[i]])[1]) text(defaults$xlim[1],kk+1.25,adj=c(0,0),paste('*Studies omitted for insufficient data:',paste(xa0[[i]],collapse=', ')),cex=.7)
    abline(h=-1,lty=2)
    addpoly(Moderator0[[i]],row=-2,mlab='SUMMARY (Random Effects)',font=2,cex=.8)
  }
  plots.mod[[i+2]] <- recordPlot()
}

pdf('Plots_Moderator.pdf',11,8.5) # start saving plots to pdf file
sapply(plots.mod,replayPlot) # generate forest plot for each analysis
graphics.off() # end saving plots

