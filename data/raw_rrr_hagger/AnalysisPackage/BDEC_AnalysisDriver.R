## Configuration: EDIT THIS SECTION BEFORE RUNNING
AnalysisDir <- '/home/kesslerd/repos/Analysis/BDEC/'
MSIT.ExperimentName <- 'MAS_MSIT_event_200_3'
LetE.Easy.ExperimentName <- 'LetterETaskt_150_3_easy'
LetE.Hard.ExperimentName <- 'LetterETaskt_150_3_hard'

### Set Working Directory (assume we're running this on freewill)
setwd(AnalysisDir)


### Functions
source('func.R')

### Libraries/dependencies
#### Install any required packages
list.of.packages <- c('reshape','car','knitr','psych','abind','zoo','retimes')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

#### Load packages
require(reshape,warn.conflicts=FALSE,quietly=TRUE)
require(car,warn.conflicts=FALSE,quietly=TRUE)
require(knitr)
require(psych)
require(abind)
require(zoo)
require(retimes)

## Loading
### Load individual files
#### Load Subject Condition Assignment Data
source('Cond_LoadParse.R')

#### Load, Parse, and Write Letter E Data
source('LetE_LoadParse.R')

#### Load, Parse, and Write MSIT Data
source('MSIT_LoadParse.R')

### Merge Files Together

SubjectData <-  subs

#### Make names more obvious
names(msit.SubjectData)[! names(msit.SubjectData) %in% names(SubjectData)] <- paste(names(msit.SubjectData)[! names(msit.SubjectData) %in% names(SubjectData)],'.MSIT',sep='')

names(LetE.SubjectData)[! names(LetE.SubjectData) %in% names(SubjectData)] <- paste(names(LetE.SubjectData)[! names(LetE.SubjectData) %in% names(SubjectData)],'.LetE',sep='')

#### Do the actual merge

SubjectData <- merge(SubjectData,LetE.SubjectData,all=TRUE,suffixes=c('','.LetE'))
SubjectData <- merge(SubjectData,msit.SubjectData,all=TRUE,suffixes=c('','.MSIT'))

## Funnel Analysis (Inclusion Criteria)

SubjectData$Include.Overall <- SubjectData$Include.Overall.MSIT & SubjectData$Include.LetE

### Recruitment N's

print('Recruitment N\'s')
print(addmargins(with(SubjectData,table(Task))))

### Usable Letter E Data
print('Included Letter E Data')
print(addmargins(with(SubjectData [ SubjectData$Include.LetE,],table(Task))))

### Usable MSIT Data
print('Included MSIT Data')
print(addmargins(with(SubjectData[ SubjectData$Include.Overall.MSIT,],table(Task))))

### Summarize Data Utilization Intersections
print('Data usability summary (Subject Cells by Inclusion Criteria')
print(addmargins(with(SubjectData,table(Task,Include.Overall.MSIT,Include.LetE)),margin=c(1,2)))

## Write Out Data

### Subset Lite version for only variables of interest

litevars <- c(
  'Subject',
  'Task',
  'Include.Overall',
  'ExGauss.I.RTVar.MSIT',
  'ExGauss.C.RTVar.MSIT')

liteSubjectData <- SubjectData [,litevars]

### Write CSVs
write.csv(SubjectData,'BDEC_Results_Full.csv',row.names=FALSE)
write.csv(liteSubjectData,'BDEC_Results_Lite.csv',row.names=FALSE)

