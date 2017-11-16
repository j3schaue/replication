## Read and Parse EDATs
### Read in file
LetE.targetfile <- dir(pattern='LetterE_EData.txt',ignore.case=T)
LetE = read.delim(LetE.targetfile,fileEncoding = 'UCS-2LE',as.is=TRUE)

### ID any missing subjects from SubjectStatus file
SubsMissing = setdiff(unique(LetE$Subject),unique(subs$Subject));
if(length(SubsMissing)!=0){
  print('I have Letter E EDATs, but no phenotypic data, for these subjects',quote=FALSE)
  print(SubsMissing)}
rm(SubsMissing)

EDATMissing = setdiff(unique(subs$Subject),unique(LetE$Subject))
if(length(EDATMissing)!=0){
  print('I have phenotypic data, but no Letter E EDATs, for these subjects',quote=FALSE)
  print(EDATMissing)}
rm(EDATMissing)

### Clean things up
#### Drop clock strings cuz they're big and ugly
LetE=LetE[,!grepl('Clock',names(LetE))]

#### Filter so that you include only tests, no practice
LetE = LetE [ LetE$ExperimentName %in% c(LetE.Easy.ExperimentName,LetE.Hard.ExperimentName),]

#### Clean up some classes
LetE$Word.RESP  = ifelse(!is.na(LetE$Word.RESP),1,0)  # recode all responses as 1, all non responses as 0
LetE$Word.CRESP = ifelse(!is.na(LetE$Word.CRESP),1,0) # recode all responses as 1, all non responses as 0

#### Recode task difficulty so strings can be shorter
LetE$TaskType = NA
LetE$TaskType[ LetE$ExperimentName==LetE.Easy.ExperimentName] = 'Easy'
LetE$TaskType[ LetE$ExperimentName==LetE.Hard.ExperimentName] = 'Hard'

#### Recode Trial Types

##### Trial Types are encoded in TrialType. 
###### 1 - No E
###### 2 - NonLonely E
###### 3 - Lonely E

##### Recode trial types also into Go/NoGo depending on TaskType
LetE$GoNoGo = NA
LetE$GoNoGo[ LetE$TaskType=='Easy' & LetE$TrialType %in% c(2,3)] = 'Go'
LetE$GoNoGo[ LetE$TaskType=='Hard' & LetE$TrialType %in% c(3)] = 'Go'
LetE$GoNoGo[ LetE$TaskType=='Easy' & LetE$TrialType %in% c(1)] = 'NoGo'
LetE$GoNoGo[ LetE$TaskType=='Hard' & LetE$TrialType %in% c(1,2)] = 'NoGo'

## Calculate Trial Response Variables
LetE=within(LetE,{   #Calculate accuracy and trial duration within the data environment (makes the code much easier to read)
  Acc= Word.RESP == Word.CRESP
  TrialDur=Word.RT
  TrialDur=TrialDur/1000
  SigDet = NA
  SigDet [Word.RESP==1 & Acc==1] = 'Hit' 
  SigDet [Word.RESP==0 & Acc==0] = 'Miss'
  SigDet [Word.RESP==1 & Acc==0] = 'False Alarm'
  SigDet [Word.RESP==0 & Acc==1] = 'Correct Rejection'
  SigDet = as.factor(SigDet)
})


## Calculate Subject-level Stats

### Confirm correct subject assignment
EDATTask = tapply(LetE$TaskType,LetE$Subject,function(x){substr(unique(x),1,1)})
EDATTask = data.frame(Subject = row.names(EDATTask),Task.EDAT = EDATTask)

diffmap = merge(subs,EDATTask)
MisMapSubs = with(diffmap,Subject [ Task != Task.EDAT])

if(length(MisMapSubs)!=0){
print('The following subjects appear to be assigned to the wrong condition',quote=FALSE)
print(MisMapSubs)
print(diffmap[ diffmap$Subject %in% MisMapSubs,c('Subject','Task','Task.EDAT')])
}

rm(
  EDATTask,
  diffmap,
  MisMapSubs
  )

### Prepare Melted Objects
#### Merge in Subject Info
LetE2=merge(LetE,subs)

#### Do the melting
LetE.melt.DurAcc=melt(LetE2,id.vars=c('Subject','Task','TrialType','GoNoGo'),measure.vars=c('TrialDur','Acc'))
LetE.melt.DurByAcc=melt(LetE2,id.vars=c('Subject','Task','TrialType','Acc','GoNoGo'),measure.vars=c('TrialDur'))
LetE.melt.SigDet = melt(LetE2,id.vars=c('Subject','Task'),measure.vars='SigDet')

### Accuracy Stats
#### Overall
LetE.OverallAcc = data.frame(cast(LetE.melt.DurAcc,Subject + Task ~ variable,fun.aggregate=mean,na.rm=TRUE,subset=variable=='Acc'))
names(LetE.OverallAcc)[3] = 'Acc.Overall'
LetE.SubjectData = LetE.OverallAcc
rm('LetE.OverallAcc')

#### Accuracy by Trial Type
LetE.AccByTrialType = data.frame(cast(LetE.melt.DurAcc,Subject + Task ~ TrialType + variable,fun.aggregate=mean,na.rm=TRUE,subset=variable=='Acc'))
names(LetE.AccByTrialType)[3:5]=c('Acc.No.E','Acc.Non.Lonely.E','Acc.Lonely.E')
LetE.SubjectData = merge(LetE.SubjectData,LetE.AccByTrialType)
rm('LetE.AccByTrialType')

#### Accuracy by Go NoGo
LetE.AccByGoNoGo = data.frame(cast(LetE.melt.DurAcc,Subject + Task ~ GoNoGo + variable,fun.aggregate=mean,na.rm=TRUE,subset=variable=='Acc'))
LetE.SubjectData = merge(LetE.SubjectData,LetE.AccByGoNoGo)
rm(LetE.AccByGoNoGo)

#### Signal Detection Stats
LetE.SigDet = data.frame(cast(LetE.melt.SigDet,Subject + Task ~ variable,fun.aggregate=table))
LetE.SubjectData = merge(LetE.SubjectData,LetE.SigDet)
rm(LetE.SigDet)

### Clean up Aggregate Stats

#### Coerce most columns to numeric
LetE.SubjectData[,3:ncol(LetE.SubjectData)] = data.frame(apply(LetE.SubjectData[,3:ncol(LetE.SubjectData)],c(1,2),function(x){as.numeric(x)}))

#### Clean up NaN Inf, etc, set it all to NA
LetE.SubjectData[,3:ncol(LetE.SubjectData)] = data.frame(apply(LetE.SubjectData[,3:ncol(LetE.SubjectData)],c(1,2),function(x){ifelse(!is.finite(x),NA,x)}))

### Subject Exclusions
LetE.SubjectData$Include <- LetE.SubjectData$Acc.Overall > .8


### Write out subject-level summary info
#write.csv(LetE.SubjectData,'LetterESummary.csv',row.names=FALSE)

## Clean up

### Remove some no longer needed variables

rm(LetE,LetE2,LetE.melt.DurAcc,LetE.melt.DurByAcc,LetE.melt.SigDet)

LetE.SubjectData$LetEPresent <- TRUE
