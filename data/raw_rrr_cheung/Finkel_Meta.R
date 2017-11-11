source('Finkel_Analysis.R') # load analysis script from source location

#====================================#
# Read and Analyze Data For Each Lab #
#====================================#

#------#
# Data #
#------#

csv <- list.files(pattern="*.csv") # names of all csv files in folder
authors <- sapply(strsplit(csv,split='_',fixed=TRUE),function(x){x[1]}) # authors' names
k <- length(authors) # number of labs
data.all <- lapply(csv,read.csv,header=FALSE) # Compile data from all labs into a single list
names(data.all) <- authors # assign author name to each dataset in list

ms <- rep(0,k)
for(l in 1:k){
  data2 <- data.all[[l]]
  colnames(data2) <- apply(data2[1,],2,as.character) # label columns
  data2 <- data2[-(1:2),] # delete first two irrelevant rows
  Exclude <- convert(data2$Exclude) # exclusion criteria
  if(max(Exclude)==2) ms[l] <- 1 # flag lab with any subjects that are Exclude=2
}
rm(Exclude,data2,l) # remove Exclude, data2, l from environment
if(max(ms)==0) m <- 1 else m <- 2 # if no labs have Exclude=2, run meta-analysis just once

for(q in 1:m){

if(q==1) include <- 0 # first run with only 0
if(q==2) include <- c(0,2) # second run with 0 and 2

#-----------------------------------------------------#
# Option 1: Run analyses and save output for each lab #
#-----------------------------------------------------#

if(q==1) {meta.temp <- vector(mode="list",length=k); ks <- 1:k} # preallocate list to store outputs and run analysis for all labs
if(q==2) ks <- which(ms==1) # run analysis for labs with Exclude=2
for(i in ks){ # for each lab
  data <- data.all[[i]] # data from lab i
  author <- authors[i] # author for lab i
  file.out <- paste(author,'_output_',paste(include,collapse=''),'.txt',sep='') # attach author to output file name
  sink(file.out) # start writing all console output to file
  cat(paste('============',paste(rep('=',nchar(author)),sep='',collapse=''),'\n',sep=''))
  cat(paste('RESULTS FOR',toupper(author)))
  cat(paste('\n============',paste(rep('=',nchar(author)),sep='',collapse=''),'\n',sep=''))
  cat(paste('\n*Inclusion:',paste(include,collapse=' & '),'\n'))
  meta.temp[[i]] <- Finkel.Analyze(data,include,10000) # analyze data and store results for meta-analysis
  sink() # end writing output
}
meta.all <- simplify2array(meta.temp) # reformat stored results to make it easier to work with


#-----------------------------------------------------------#
# Option 2: Run analyses without saving output for each lab #
#-----------------------------------------------------------#

#meta.all <- sapply(data.all,Finkel.Analyze,include,100) # analyze data from each lab and store results for meta-analysis


#===============#
# Meta-Analyses #
#===============#

if(suppressWarnings(!require(metafor,quietly=TRUE))) install.packages('metafor') # install metafor package if necessary
library(metafor) # load metafor package

for(i in 1:nrow(meta.all)){ # for each dependent variable
  
  #-------#
  # Setup #
  #-------#
  
  # Data #
  label.DV <- row.names(meta.all)[i] # label for current dependent variable
  if(label.DV=='SubjCommit') name.DV <- 'Subjective Commitment' else name.DV <- paste(label.DV,'Forgiveness')
  nr <- nrow(meta.all[[i]]) # number of analyses
  DV <- array(unlist(meta.all[i,]),c(nr,6,k))
    # arrange results for dependent variable from all labs into 3D array
    # dimension 1: type of analysis (t, ANOVA, ANCOVA, etc.)
    # dimension 2: results (mean for high, mean for low, difference, SE)
    # dimension 3: lab
  dimnames(DV) <- c(dimnames(as.data.frame(meta.all[[i]])),list(authors)) # label dimensions
  
  # RMA Output #
  rma.out <- paste(label.DV,'_rma_',paste(include,collapse=''),'.txt',sep='') # attach author to output file name
  sink(rma.out) # start writing rma output for current DV
  cat(paste('\n========',paste(rep('=',nchar(label.DV)),sep='',collapse=''),'\n',sep=''))
  cat(paste('RMA FOR',toupper(label.DV)))
  cat(paste('\n========',paste(rep('=',nchar(label.DV)),sep='',collapse=''),'\n',sep=''))
  cat(paste('\n*Inclusion:',paste(include,collapse=' & '),'\n\n'))
  
  # Forest Plot Output #
  forest.plots <- vector('list',nr) # allocate list for storing plots
  par(mar=c(6,6,3,4)) # set margins for plot window
  
  #----------#
  # Analyses #
  #----------#
  
  DV.rma <- vector('list',nr)
  analysis.type <- c('Main','Overall','Male','Female','Overall','Male','Female')
  mod.med <- c('Moderator: None',rep('Moderator: Gender',6))
  covar <- c(rep('Covariates: None',4),rep('Covariates: Impression Management (IM), Self-Deception (SD)',3))
  for(j in 1:nr){ # for each analysis
    
    # Data #
    H <- round(DV[j,1,],2) # means for high across labs (rounded to 2 decimal places)
    L <- round(DV[j,2,],2) # means for low across labs (rounded to 2 decimal places)
    Diff <- DV[j,3,] # effect size across labs
    SE <- DV[j,4,] # standard errors across labs
    K <- DV[j,5,] # number of items across labs
    N <- DV[j,6,]
    
    # RMA #
    cat(paste(rep('-',70,sep='',collapse='')),'\n\n',sep='')
    if(j<8) cat('[',analysis.type[j],': High-Low]','\n\n',mod.med[j],'\n',covar[j],'\n',sep='')  
    if(j==8 & nr==9) cat('[Mediation: a*b]\n\n','Mediator: Subjective Commitment\n',sep='')
    if((j==8 & nr==8) | (j==9 & nr==9)) cat("[Reliability: Cronbach's Alpha]\n\n",'Missing Data: Listwise Deletion\n',sep='')
    if(j<8 | (j==8 & nr==9)) DV.rma[[j]] <- rma(yi=Diff,sei=SE,slab=authors,knha=TRUE) 
    if((j==8 & nr==8) | (j==9 & nr==9)) DV.rma[[j]] <- rma(ai=Diff,mi=K,ni=N,measure='ARAW',slab=authors,knha=TRUE)
      # random effects meta-analyses with Knapp & Hartung Adjustment
    print(DV.rma[[j]],signif.legend=FALSE,signif.stars=FALSE) # print results without significance codes
    
    # Forest Plot (Individual) #
    if(j>1){
      defaults <- forest(DV.rma[[j]]) # default forest plot
      if(j<8){
        xp1 <- defaults$xlim[1]+.175*(defaults$xlim[2]-defaults$xlim[1]) # x-coordinate for high means column
        xp2 <- defaults$xlim[1]+.255*(defaults$xlim[2]-defaults$xlim[1]) # x-coordinate for low means column
        DV.forest <- forest(DV.rma[[j]],ilab=cbind(H,L),ilab.xpos=c(xp1,xp2),xlab='Difference in Means',addfit=FALSE,ylim=c(-2,k+5))
          # forest plot of results
          # meta-analysis summary omitted (edited version to be added below)
          # means for high and low included for each lab
        text(xp1,k+1.08,'High',font=2) # column name for high means
        text(xp2,k+1.11,'Low',font=2) # column name for low means
        text(defaults$xlim[2],k+1,adj=c(1.05,0),'High - Low [95% CI]',font=2) # column name for confidence intervals
        text(mean(defaults$xlim),k+5,paste(name.DV,' [',analysis.type[j],']',sep=''),font=2,cex=1.25) # title indicating variable of interest
        text(mean(defaults$xlim),k+4,mod.med[j],font=2,cex=1) # subtitle indicating moderator
        text(mean(defaults$xlim),k+3.3,covar[j],font=2,cex=1) # subtitle indicating covariates
      }
      if(j==8 & nr==9){
        DV.forest <- forest(DV.rma[[j]],xlab='Mediation Effect',addfit=FALSE,ylim=c(-2,k+5))
          # forest plot of results
          # meta-analysis summary omitted (edited version to be added below)
          # means for high and low included for each lab
        text(defaults$xlim[2],k+1,adj=c(1.1,0),'a*b [95% CI]',font=2) # column name for confidence intervals
        text(mean(defaults$xlim),k+4.5,paste(name.DV,' [Mediation]',sep=''),font=2,cex=1.25) # title indicating variable of interest
        text(mean(defaults$xlim),k+3.5,'Mediator: Subjective Commitment',font=2,cex=1) # subtitle indicating mediator
        text(mean(defaults$xlim),k+3.3,'',font=2,cex=1) # subtitle indicating covariates
      }
      if((j==8 & nr==8) | (j==9 & nr==9)){
        DV.forest <- forest(DV.rma[[j]],xlab="Cronbach's Alpha",addfit=FALSE,ylim=c(-2,k+5))
        # forest plot of results
        # meta-analysis summary omitted (edited version to be added below)
        # means for high and low included for each lab
        text(defaults$xlim[2],k+1,adj=c(1.1,0),expression(bold(paste(alpha,' [95% CI]')))) # column name for confidence intervals
        text(mean(defaults$xlim),k+4,paste(name.DV,' [Reliability]',sep=''),font=2,cex=1.25) # title indicating variable of interest
        text(defaults$xlim[2],k+4,adj=c(1.05,0),'(Listwise Deletion)',cex=.75) # inclusion criteria
        text(mean(defaults$xlim),k+3.3,'',font=2,cex=1) # subtitle indicating covariates
      }
      addpoly(DV.rma[[j]],row=-1,mlab='SUMMARY (Random Effects)',font=2) 
        # replace meta-analysis summary with this version (remaned and bolded)
      abline(h=0) # horizontal line to separate meta-analysis summary from individual results
      abline(h=k+3,col='white')
      segments(0,k+2.01,0,k+3,col='white')
      text(defaults$xlim[1],k+1,adj=c(-0.15,0),'Study',font=2) # column name for labs
      text(defaults$xlim[2],k+5,adj=c(1.1,0),paste('*Inclusion:',paste(include,collapse=' & ')),cex=.75)
      abline(h=k+2)
      forest.plots[[j]] <- recordPlot()
    }
    
  }
  
  # Forest Plot (Combined) #
  defaults <- forest(DV.rma[[1]]) 
  xp1 <- defaults$xlim[1]+.175*(defaults$xlim[2]-defaults$xlim[1])
  xp2 <- defaults$xlim[1]+.255*(defaults$xlim[2]-defaults$xlim[1])
  H <- round(DV[1,1,],2)
  L <- round(DV[1,2,],2) 
  DV.forest <- forest(DV.rma[[1]],ilab=cbind(H,L),ilab.xpos=c(xp1,xp2),xlab='Difference in Means',addfit=FALSE,ylim=c(-13,k+5),cex=.75)
  abline(h=k+3,col='white')
  segments(0,k+2.01,0,k+3,col='white')
  text(xp1,k+1.15,'High',font=2,cex=.75) 
  text(xp2,k+1.21,'Low',font=2,cex=.75)
  text(defaults$xlim[2],k+1,adj=c(1.05,0),'High - Low [95% CI]',font=2,cex=.75) 
  text(defaults$xlim[1],k+1,adj=c(-0.15,0),'Study',font=2,cex=.75)
  text(mean(defaults$xlim),k+5,paste(name.DV,' [',analysis.type[1],']',sep=''),font=2,cex=1.25) 
  text(mean(defaults$xlim),k+3.9,mod.med[1],font=2,cex=1)
  text(mean(defaults$xlim),k+2.9,covar[1],font=2,cex=1) 
  abline(h=k+2)
  text(defaults$xlim[2],k+5,adj=c(1.1,0),paste('*Inclusion:',paste(include,collapse=' & ')),cex=.75)
  abline(h=0) 
  text(defaults$xlim[1],-1,adj=c(0,0),'SUMMARY (Random Effects)',font=2,cex=.75) 
  addpoly(DV.rma[[1]],row=-2,mlab='Main',font=2,cex=.75)
  abline(h=-3,lty=2)
  text(defaults$xlim[1],-4,adj=c(-.07,0),'Gender Moderator',font=2,cex=.75)
  mlabs <- c('Overall','Male','Female')
  for(j in 2:4) addpoly(DV.rma[[j]],row=-(j+3),mlab=mlabs[j-1],cex=.75)
  abline(h=-8,lty=2)
  text(defaults$xlim[1],-9,adj=c(-.03,0),'Gender Moderator, IM & SD Covariates',font=2,cex=.75)
  for(j in 5:7) addpoly(DV.rma[[j]],row=-(j+5),mlab=mlabs[j-4],cex=.75)
  forest.plots[[1]] <- recordPlot()
  
  # Output Files #
  sink() # end writing rma output for current DV
  forest.out <- paste(label.DV,'_forest_',paste(include,collapse=''),'.pdf',sep='') # attach DV label to output file name
  pdf(forest.out,11,8.5) # start saving plots to pdf file
  sapply(forest.plots,replayPlot) # generate forest plot for each analysis
  graphics.off() # end saving plots
  
}

}

#==================#
# Demographic Info #
#==================#

demo <- t(sapply(data.all,Finkel.Demo)) # run demographic analysis for all labs
colnames(demo) <- c('N_Total','N_Male','N_Female','Age_Mean','Age_SD','Male_Exclude_1','Female_Exclude_1','Male_Exclude_1&2','Female_Exclude_1&2')
write.csv(demo,'Demographics.csv')

#==============#
# Correlations #
#==============#

cortabs <- vector('list',2)
metacor <- matrix(0,2,6)

for(i in 1:2){ # for each exclusion rule
  
  if(i==1) include=0 else include=c(0,2) # inclusion criteria
  
  # Correlations #
  cortab <- t(sapply(data.all,Finkel.Cor,include)) # run correlation analysis for each lab
  pair.DV <- c('Exit-Neglect','Exit-Voice','Exit-Loyalty','Neglect-Voice','Neglect-Loyalty','Voice-Loyalty','N_Complete') # column names 
  colnames(cortab) <- pair.DV # label columns
  cortabs[[i]] <- cortab
  
  # Setup #
  rma.out <- paste('Correlations_rma_',paste(include,collapse=''),'.txt',sep='') # attach DV pair to output file name
  sink(rma.out) # start writing rma output for current DV
  cat(paste('\n====================\n'))
  cat(paste('RMA FOR CORRELATIONS'))
  cat(paste('\n====================\n'))
  cat(paste("\nResults in Fisher's z\n(Listwise Deletion)\n"))
  cat(paste('\n*Inclusion:',paste(include,collapse=' & '),'\n\n'))
  forest.plots <- vector('list',6) # allocate list for storing plots
  par(mar=c(6,6,3,4)) # set margins for plot window
  
  for(j in 1:6){
  
  # RMA #
  cat(paste(rep('-',70,sep='',collapse='')),'\n\n',sep='')
  cat('[',pair.DV[j],']','\n',sep='')  
  rma.cor <- rma(ri=cortab[,j],ni=cortab[,7],measure='ZCOR',slab=authors,knha=TRUE) 
    # random effects meta-analysis with Knapp & Hartung Adjustment
  print(rma.cor,signif.legend=FALSE,signif.stars=FALSE) # print results without significance codes
  metacor[i,j] <- round(predict(rma.cor,transf=transf.ztor)$pred,2) # summary correlation
  
  # Forest Plot #
  defaults <- forest(rma.cor,transf=transf.ztor) # default forest plot
  forest.cor <- forest(rma.cor,transf=transf.ztor,addfit=FALSE,ylim=c(-2,k+4))
    # forest plot of results
    # meta-analysis summary omitted (edited version to be added below)
  text(defaults$xlim[2],k+1,adj=c(1.1,0),'r [95% CI]',font=2) # column name for confidence intervals
  text(mean(defaults$xlim),k+3.5,paste(pair.DV[j],'Correlation'),font=2,cex=1.25) # title indicating variable of interest
  addpoly(rma.cor,transf=transf.ztor,row=-1,mlab='SUMMARY (Random Effects)',font=2) 
    # replace meta-analysis summary with this version (remaned and bolded)
  abline(h=0) # horizontal line to separate meta-analysis summary from individual results
  text(defaults$xlim[1],k+1,adj=c(-0.15,0),'Study',font=2) # column name for labs
  text(defaults$xlim[2],k+3.5,adj=c(1.1,0),paste('*Inclusion:',paste(include,collapse=' & ')),cex=.75) # inclusion criteria
  text(defaults$xlim[2],k+3,adj=c(1.05,0),'(Listwise Deletion)',cex=.75) # inclusion criteria
  forest.plots[[j]] <- recordPlot()
  
  }
  
  # Output Files #
  sink() # end writing rma output for current DV
  forest.out <- paste('Correlations_forest_',paste(include,collapse=''),'.pdf',sep='') # attach DV label to output file name
  pdf(forest.out,11,8.5) # start saving plots to pdf file
  sapply(forest.plots,replayPlot) # generate forest plot for each analysis
  graphics.off() # end saving plots

}

for(i in 1:2){

cor.l <- matrix('',4,4)
cor.l[lower.tri(cor.l)] <- metacor[i,]
cor.l <- data.frame(cor.l[,-4])
dimnames(cor.l) <- list(c('Exit','Neglect','Voice','Loyalty'),c('Exit','Neglect','Voice'))
if(i==1) {
  sink('Correlation_Matrix_0.txt') 
  cat('\n*Inclusion: 0\n')
} else {
  sink('Correlation_Matrix_02.txt')
  cat('\n*Inclusion: 0 & 2\n')
}
cat
cat('\n(Listwise Deletion)\n\n')
cat('\n---------------------------')
cat('\nIndividual Lab Correlations')
cat('\n---------------------------\n\n')
print(cortabs[[i]])
cat('\n\n-----------------------')
cat('\nMeta Correlation Matrix')
cat('\n-----------------------\n\n')
print(cor.l)
sink()

}