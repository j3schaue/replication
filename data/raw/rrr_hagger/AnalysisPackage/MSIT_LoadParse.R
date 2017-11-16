## Read and Parse EDATs
### Read in file
msit.targetfile <- dir(pattern='MSIT200_EData.txt',ignore.case=T)
msit = read.delim(msit.targetfile,fileEncoding = 'UCS-2LE',as.is=TRUE)

### ID any missing subjects from SubjectStatus file
SubsMissing = setdiff(unique(msit$Subject),unique(subs$Subject))
if(length(SubsMissing)!=0){
  print('I have MSIT EDATs, but no phenotypic data, for these subjects',quote=FALSE)
  print(SubsMissing)}
rm(SubsMissing)

EDATMissing = setdiff(unique(subs$Subject),unique(msit$Subject))
if(length(EDATMissing)!=0){
  print('I have phenotypic data, but no MSIT EDATs, for these subjects',quote=FALSE)
  print(EDATMissing)}
rm(EDATMissing)

### Clean things up
#### Drop clock strings cuz they're big and ugly
msit=msit[,!grepl('Clock',names(msit))]

#### Filter so that you include only tests, no practice
msit=msit[msit$ExperimentName==MSIT.ExperimentName,]

#### Clean up some classes
msit$Digits.RESP <- as.numeric(msit$Digits.RESP) # coerce responses to be numeric
msit$Digits.RESP[is.na(msit$Digits.RESP)] <- 0 # set non responses to zero
msit$Digits.CRESP <- msit$target # nonresponse/miss trials result in CRESP being NA, so overwrite with target
msit$Digits.CRESP <- as.numeric(msit$Digits.CRESP) # ensure it is numeric


## Calculate trial response variables
msit=within(msit,{   #Calculate accuracy and trial duration within the data environment (makes the code much easier to read)
  Acc = as.numeric(Digits.RESP == Digits.CRESP)
  TrialDur=Digits.RT
  TrialDur=TrialDur/1000
  TrialTypeNum=TrialType
  TrialTypeNumAccOnly=ifelse(!is.na(TrialTypeNum) & is.na(TrialDur),max(TrialTypeNum,na.rm=TRUE)+2,TrialTypeNum) # recode omission errors to 4
  TrialTypeNumAccOnly=ifelse(!is.na(TrialTypeNum) & !is.na(TrialDur) & Acc==0,max(TrialTypeNum,na.rm=TRUE)+1,TrialTypeNumAccOnly) # recode comission errors to 3
})

## Compute Subject-Level Stats

### Prepare the melted objects

#### Merge in Phenotypic Data
msit2=merge(msit,subs)

#### Do the melting
msit.melt=melt(msit2,id.vars=c('Subject','Task','TrialTypeNum','Acc'),measure.vars='TrialDur')
msit.accmelt <- melt(msit2,id.vars=c('Subject','Task','TrialTypeNum'),measure.vars='Acc')
##### Recode Trial Types to be more user friendly
msit.melt$TrialTypeNum[msit.melt$TrialTypeNum==1]='C'
msit.melt$TrialTypeNum[msit.melt$TrialTypeNum==2]='I'
msit.accmelt$TrialTypeNum[msit.accmelt$TrialTypeNum==1]='C'
msit.accmelt$TrialTypeNum[msit.accmelt$TrialTypeNum==2]='I'


### Do some casting

#### RT Metrics
##### Means by:
###### Trial Type (Accurate Only)
msit.cast.means=data.frame(cast(msit.melt,Subject + Task  ~  TrialTypeNum + Acc,fun.aggregate=mean,subset=(Acc==1)))
names(msit.cast.means)[3:4] = c('C_1_MeanRT','I_1_MeanRT')
msit.SubjectData = msit.cast.means

###### Accurate Trials (Across Trial Type)
msit.cast.means.overall=data.frame(cast(msit.melt,Subject + Task  ~   Acc,fun.aggregate=mean,subset=(Acc==1)))
names(msit.cast.means.overall)[3]='Overall_1_RT'
msit.SubjectData <- merge(msit.SubjectData,msit.cast.means.overall)
rm(msit.cast.means.overall)
##### RT SD by:
###### Trial Type (Accurate Only)
msit.cast.sd=data.frame(cast(msit.melt,Subject + Task  ~  TrialTypeNum + Acc,fun.aggregate=sd,subset=(Acc==1)))
names(msit.cast.sd)[3:4] = c('C_1_SDRT','I_1_SDRT')
msit.SubjectData = merge(msit.SubjectData,msit.cast.sd)
rm(msit.cast.sd)

###### Accurate Trials (Across Trial Type)
msit.cast.sd.overall=data.frame(cast(msit.melt,Subject + Task  ~  Acc,fun.aggregate=sd,subset=(Acc==1)))
names(msit.cast.sd.overall)[3]='Overall_1_SD'
msit.SubjectData <- merge(msit.SubjectData,msit.cast.sd.overall)
rm(msit.cast.sd.overall)

#### Accuracy Metrics
##### by Trial Type
msit.cast.acc = data.frame(cast(msit.accmelt,Subject + Task ~ TrialTypeNum, fun.aggregate = mean))
names(msit.cast.acc)[3:4] <- paste('Acc.',names(msit.cast.acc)[3:4],sep='')

msit.SubjectData <- merge(msit.SubjectData,msit.cast.acc)
rm(msit.cast.acc)

##### Across Trial Types (Overall)
msit.cast.acc = data.frame(cast(msit.accmelt,Subject + Task ~ ., fun.aggregate = mean))
names(msit.cast.acc)[3] = 'Acc.Overall'
msit.SubjectData <- merge(msit.SubjectData,msit.cast.acc)
rm(msit.cast.acc)

#### Reponse vs Nonresponse Metrics
msit.cast.response = data.frame(cast(msit.melt,Subject + Task ~ TrialTypeNum, fun.aggregate= function(x) mean(x==0)))
names(msit.cast.response)[3:4] = c('C_NonResponse','I_NonResponse')
msit.SubjectData <- merge(msit.SubjectData,msit.cast.response)
rm(msit.cast.response)

### Calculate Derived Subjectlevel Stats
#### Trial Type Effect (I minus C)
##### RT Mean
msit.SubjectData$TrialTypeEffect_Mean = msit.SubjectData$I_1_MeanRT - msit.SubjectData$C_1_MeanRT

##### RT SD
msit.SubjectData$TrialTypeEffect_SD = msit.SubjectData$I_1_SDRT - msit.SubjectData$C_1_SDRT


## RT Cleaning Pre Fitting

### Melt dataset down
names(msit2)[names(msit2)=='TrialNum'] = 'TrialNumber'

msit.spec.melt=melt(msit2,id.vars=c('Subject','Task','TrialNumber','TrialTypeNum'),measure.vars=c('TrialDur','Acc','TrialType'))

msit.spec.cast = cast(msit.spec.melt,Subject ~ TrialNumber ~ variable)

### Preprocess RT's

#### Censor inaccurate trials

msit.spec.cast = abind(msit.spec.cast, ifelse(msit.spec.cast[,,2]==1,msit.spec.cast[,,1],NA),along=3)
dimnames(msit.spec.cast)[[3]][4] = 'TrialDur.Censored'


## RT Distribution Fitting
### Ex Gaussian Fit
#### Get the Censored Trial Durations (Generated up above for spectral estimates)
CensoredTrialDur = msit.spec.cast[,,c('TrialDur.Censored','TrialType')]

#### Estimate Mu, Sigma, and Tau for each subject (for available RTs)
ExGauss <- aperm(apply(CensoredTrialDur,1,function(x){
  x <- x[!is.na(x[,1]),]
  Type1 <- x[ x[,2]==1,1]
  Type2 <- x[ x[,2]==2,1]
  if(length(Type1)>2){
    Type1Ex <- mexgauss(Type1)
  } else {
    Type1Ex <- c(mu=NA,sigma=NA,tau=NA)
  }
  if(length(Type2)>2) {
    Type2Ex <- mexgauss(Type2)
  } else{
    Type2Ex <- c(mu=NA,sigma=NA,tau=NA)
  }
  ExG <- c(Type1Ex,Type2Ex)
  names(ExG) <- c('C.mu','C.sigma','C.tau','I.mu','I.sigma','I.tau')
  return(ExG)
}))
ExGauss <- data.frame(Subject=rownames(ExGauss),ExGauss=ExGauss)

ExGauss$ExGauss.C.RTVar <- ExGauss$ExGauss.C.sigma + ExGauss$ExGauss.C.tau
ExGauss$ExGauss.I.RTVar <- ExGauss$ExGauss.I.sigma + ExGauss$ExGauss.I.tau

#### Bind it back onto msit.SubjectData
msit.SubjectData <- merge(msit.SubjectData,ExGauss,all=TRUE)


## Subject Exclusions

### Number of Errors (greater than 80% accuracy)
msit.SubjectData$Include.Errors <- msit.SubjectData$Acc.Overall > .8

### 2SD Filter

#### Calc Grand Mean and SD for Incongruent Correct
I.1.MeanRT = mean(msit.SubjectData$I_1_MeanRT,na.rm=TRUE)
I.1.SD = sd(msit.SubjectData$I_1_MeanRT,na.rm=TRUE)

#### Calc Grand Mean and SD for Incongruent RTV
I.1.RTV.Mean <- mean(msit.SubjectData$ExGauss.I.RTVar,na.rm=TRUE)
I.1.RTV.SD <-  sd(msit.SubjectData$ExGauss.I.RTVar,na.rm=TRUE)


### Check for Subjects outside this
msit.SubjectData$Include.I2SD = msit.SubjectData$I_1_MeanRT <= I.1.MeanRT + 2*I.1.SD
msit.SubjectData$Include.IRTV2SD = msit.SubjectData$ExGauss.I.RTVar <= I.1.RTV.Mean + 2*I.1.RTV.SD


rm(
  I.1.MeanRT,
  I.1.SD,
  I.1.RTV.Mean,
  I.1.RTV.SD
  )

### Complete Cell Info
msit.SubjectData$Include.Task = msit.SubjectData$Task %in% c('E','H')

### Overal Inclusion

msit.SubjectData$Include.Overall = with(msit.SubjectData,Include.Errors & Include.I2SD & Include.IRTV2SD & Include.Task )


## Clean up the workspace
rm(
  msit,
  msit2,
  msit.cast.means,
  msit.melt,
  msit.spec.melt,
  msit.spec.cast
  )

msit.SubjectData$Include.MSITPresent <- TRUE

### Write out the file
#write.csv(msit.SubjectData,'MSITSummary.csv',row.names=FALSE)
