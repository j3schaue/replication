## N.B.: The scales used for producing the figures may change
## for the analysis of the real data, however, the core
## analysis will not change. 

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../raw/rrr_wagenmakers/")

##@@ CLEAR WORKSPACE @@##

rm(list = ls())


##@@ REQUIRED PACKAGES @@##

# function to load/install packages (credit: CousinCocaine on Stack Exchange)
packages <- function(x) {
  
  x <- as.character(match.call()[[2]])
  
  if ( ! require(x, character.only = TRUE)) {
    
    install.packages(pkgs = x, repos = "http://cran.r-project.org")
    require(x, character.only = TRUE)    
  }
}

packages(BayesFactor)
packages(compute.es)
packages(metafor)
packages(MCMCpack)
library(dplyr)

##@@ load Josine Verhagen's Replication Bayes Factor Functions @@##

source("Analysis_R_Scripts//Repfunctionspack.R")

##@@ function to convert column to numeric @@##

convert2num <- function(column) {
  
  if (is.character(column)) {
    
    column.new <- ifelse(column == "", NA, column)
    as.numeric(column.new)
    
  } else if (is.factor(column)) {
    
    as.numeric(as.character(column)) # should not be needed since "stringsAsFactors = FALSE"
    
  } else if (is.numeric(column)) {
    
    column
    
  }
}

facialFeedbackAnalysis <- function(filename, excludeParticipants) {
  #
  #  function that analyzes a .csv file with data of a facial feedback hypothesis replication experiment that is in the format of "TemplateDatafile.xlsx"
  #  and has been converted to .csv (with added second column: Participant ID)
  #
  #  input: filename = filename of the .csv data file
  #         excludeParticipants: logical -- indicates whether participants that do not meet the criteria for
  #                              inclusion are excluded
  #
  #  output: named list with the following elements:
  #
  #         ratingsSmile: participants' cartoon ratings in the smile condition
  #         ratingsPout: participants' cartoon ratings in the pout condition
  #         nSmile: number of participants in the smile condition
  #         nSmileExcluded: number of participants in the smile condition that were excluded
  #         nPout: number of participants in the pout condition
  #         nPoutExcluded: number of participants in the pout condition that were excluded
  #         meanSmile: mean cartoon rating in the smile condition
  #         meanPout: mean cartoon rating in the pout condition
  #         meanRating: overall mean rating of participants (i.e., based on participants in both conditions)
  #         sdSmile: standard deviation of the cartoon ratings in the smile condition
  #         sdPout: standard deviation of the cartoon ratings in the pout condition
  #         tValue: t-value of a one-sided independent samples t-test (equal variances assumed) which tests the hypothesis that the mean rating in the smile condition is larger
  #         pValue: p-value of the t-test
  #         df: degrees of freedom of the t-test
  #         BFplus0: independent samples t-test Bayes factor in favor of the hypothesis that the smile condition has a larger mean rating ( cauchy prior width = 1/sqrt(2) )
  #         BFr0: replication Bayes factor as proposed in Verhagen & Wagenmakers (2014)
  #         d: Cohen's d effect size
  #         dSE: standard error of d
  #         dLowerCI: lower bound of a 95% confidence interval for d
  #         dUpperCI: upper bound of a 95% confidence interval for d
  #         g: Hedges' g effect size
  #         gSE: standard error of g
  #         gLowerCI: lower bound of a 95% confidence interval for g
  #         gUpperCI: upper bound of a 95% confidence interval for g
  #         rawMeanDiff: difference between the mean of the smile and the pout condition
  #         rawMeanDiffSE: standard error of difference between the mean of smile and pout condition
  #         dataRaw: raw data without any manipulation
  #         dataAnalysis: data that were used for the analysis
  #         filename: name of the data file
  #
  
  ##@@ DATA LOADING @@##
  dataRawTmp <- read.csv(file = filename, header = TRUE, stringsAsFactors = FALSE)
  dataRaw <- dataRawTmp[-1, ] # remove second header row
  
  ##@@ DATA MANIPULATION @@##
  
  # rename columns
  colnames(dataRaw) <- c("subjectNo", "participantID", "condition", "performedCorrectlyCartoon1",
                         "performedCorrectlyCartoon2", "performedCorrectlyCartoon3",
                         "performedCorrectlyCartoon4", "performedCorrectlyTotal",
                         "ratingTask1", "ratingTask2", "ratingCartoon1",
                         "ratingCartoon2", "ratingCartoon3", "ratingCartoon4",
                         "selfReportedPerformance", "comprehensionCartoons",
                         "awareOfGoal", "participantsGuessedGoal", "age", "gender",
                         "student", "occupationFieldOfStudy")
  
  # make sure that the correct columns are numeric
  dataRaw$condition <- convert2num(dataRaw$condition)
  dataRaw$performedCorrectlyTotal <- convert2num(dataRaw$performedCorrectlyTotal)
  dataRaw$performedCorrectlyCartoon1 <- convert2num(dataRaw$performedCorrectlyCartoon1)
  dataRaw$performedCorrectlyCartoon2 <- convert2num(dataRaw$performedCorrectlyCartoon2)
  dataRaw$performedCorrectlyCartoon3 <- convert2num(dataRaw$performedCorrectlyCartoon3)
  dataRaw$performedCorrectlyCartoon4 <- convert2num(dataRaw$performedCorrectlyCartoon4)
  dataRaw$ratingCartoon1 <- convert2num(dataRaw$ratingCartoon1)
  dataRaw$ratingCartoon2 <- convert2num(dataRaw$ratingCartoon2)
  dataRaw$ratingCartoon3 <- convert2num(dataRaw$ratingCartoon3)
  dataRaw$ratingCartoon4 <- convert2num(dataRaw$ratingCartoon4)
  dataRaw$comprehensionCartoons <- convert2num(dataRaw$comprehensionCartoons)
  dataRaw$awareOfGoal <- convert2num(dataRaw$awareOfGoal)
  
  # remove rows without cartoon ratings
  cartoonRatingSubset1 <- dataRaw[ , c("ratingCartoon1", "ratingCartoon2", "ratingCartoon3", "ratingCartoon4")]
  noObsIndex <- apply(cartoonRatingSubset1, 1, function(row) all(is.na(row)) )
  dCleaned1 <-  dataRaw[ ! noObsIndex, ]
  
  nSmileBeforeExclusion <- sum(dCleaned1$condition == 1)
  nPoutBeforeExclusion <- sum(dCleaned1$condition == 0)
  
  if (excludeParticipants) {
    
    ##@@ CHECK EXCLUSION CRITERIA: 1 @@##
    
    # check whether participants performed 3 or 4 cartoon tasks correctly
    performedCorrectlyIndex <- ! is.na(dCleaned1$performedCorrectlyTotal) & dCleaned1$performedCorrectlyTotal >= 3    
    dCleaned2 <- dCleaned1[performedCorrectlyIndex, ]
    
    # check whether participants understood the cartoons
    comprehensionCartoonsIndex <- ! is.na(dCleaned2$comprehensionCartoons) &  dCleaned2$comprehensionCartoons == 1
    dCleaned3 <- dCleaned2[comprehensionCartoonsIndex, ]
    
    # check whether participants were aware of goal
    notAwareOfGoalIndex <- ! is.na(dCleaned3$awareOfGoal) & dCleaned3$awareOfGoal == 0
    dCleaned4 <- dCleaned3[notAwareOfGoalIndex, ]
    
    performedCorrectlyIndices <- dCleaned4[ , c("performedCorrectlyCartoon1", "performedCorrectlyCartoon2",
                                                "performedCorrectlyCartoon3", "performedCorrectlyCartoon4")]
    
  } else {
    
    dCleaned4 <- dCleaned1
    performedCorrectlyIndices <- matrix(1, nrow = nrow(dCleaned4), ncol = 4)
  }
  
  # calculate mean rating for each participant
  cartoonRatingSubset2 <- dCleaned4[ , c("ratingCartoon1", "ratingCartoon2", "ratingCartoon3", "ratingCartoon4")]
  
  meanCartoonRating <- numeric( nrow(cartoonRatingSubset2) )
  
  for (i in seq_len( nrow(cartoonRatingSubset2) ))
    meanCartoonRating[i] <- mean( as.numeric(cartoonRatingSubset2[i, which(performedCorrectlyIndices[i, ] == 1)]), na.rm = TRUE )
  
  dCleaned4$meanCartoonRating <- meanCartoonRating
  
  ratingsSmile <- dCleaned4[dCleaned4$condition == 1, "meanCartoonRating"]
  ratingsPout <- dCleaned4[dCleaned4$condition == 0, "meanCartoonRating"]
  
  if (excludeParticipants) {
    
    ##@@ CHECK EXCLUSION CRITERIA: 2 @@##
    
    # exclude participants that are more than 2.5 standard deviations away from condition mean
    
    outliersSmile <- (dCleaned4$condition == 1) & ((dCleaned4$meanCartoonRating > (mean(ratingsSmile) + 2.5 * sd(ratingsSmile)))
                                                   | (dCleaned4$meanCartoonRating < (mean(ratingsSmile) - 2.5 * sd(ratingsSmile))))
    
    outliersPout <- (dCleaned4$condition == 0) & ((dCleaned4$meanCartoonRating > (mean(ratingsPout) + 2.5 * sd(ratingsPout)))
                                                  | (dCleaned4$meanCartoonRating < (mean(ratingsPout) - 2.5 * sd(ratingsPout))))
    
    dCleaned4 <- dCleaned4[ ! outliersSmile & ! outliersPout, ]
    
    ratingsSmile <- dCleaned4[dCleaned4$condition == 1, "meanCartoonRating"]
    ratingsPout <- dCleaned4[dCleaned4$condition == 0, "meanCartoonRating"]
    
  }
  
  nSmile <- length(ratingsSmile)
  nSmileExcluded <- nSmileBeforeExclusion - nSmile
  meanSmile <- mean(ratingsSmile)
  sdSmile <- sd(ratingsSmile)
  
  nPout <- length(ratingsPout)
  nPoutExcluded <- nPoutBeforeExclusion - nPout
  meanPout <- mean(ratingsPout)
  sdPout <- sd(ratingsPout)
  
  # overall meanRating
  meanRating <- mean(c(ratingsSmile, ratingsPout))
  
  #@@ DATA ANALYSIS @@##
  
  # classical independent samples t-test (one-tailed, equal variances assumed)
  ttestClassical <- t.test(ratingsSmile, ratingsPout, alternative = "greater", var.equal = TRUE)
  tValue <- unname(ttestClassical$statistic)
  pValue <- ttestClassical$p.value
  df <- unname(ttestClassical$parameter)
  
  # Bayesian independent samples t-test (one-tailed)
  ttestBayesian <- ttestBF(ratingsSmile, ratingsPout, nullInterval = c(0, Inf))
  BFplus0 <- extractBF(ttestBayesian, onlybf=TRUE)[1]
  
  # Replication Bayes factor as proposed in Verhagen & Wagenmakers (2014)
  #
  # since the original article does not report relevant independent samples t-test nor
  # the number of participants, we assume equal number of participants in smile and pout
  # condition, and equal standard deviations in all three conditions of original experiment
  # (smile, neutral, pout)
  n.smile.original <- 32
  n.pout.original <- 32
  t.original <- 1.85
  t.replication <- tValue
  BFr0 <- ReplicationBayesfactorNCT(t.original, t.replication, n1 = n.smile.original, n2 = nSmile,
                                    m1 = n.pout.original, m2 = nPout, sample = 2)$BF
  BFr0 <- unname(BFr0) # remove "t" label
  
  ### calculate effect sizes and confidence intervals ###
  
  silent <- capture.output( effectSize <- tes(t=tValue, n.1=nSmile, n.2=nPout) ) # capture.output allows to suppress printing to console
  
  # Cohen's d
  d <- effectSize$d
  dSE <- sqrt(effectSize$var.d)
  dLowerCI <- effectSize$l.d
  dUpperCI <- effectSize$u.d
  
  # Hedges' g
  g <- effectSize$g
  gSE <- sqrt(effectSize$var.g)
  gLowerCI <- effectSize$l.g
  gUpperCI <- effectSize$u.g
  
  # raw mean difference
  rawMeanDiff <- meanSmile - meanPout
  rawMeanDiffSE <- rawMeanDiff/tValue
  
  #@@ STORE DATA AND ANALYSIS @@##
  studyResult <- list(ratingsSmile=ratingsSmile, ratingsPout=ratingsPout, nSmile=nSmile, nSmileExcluded=nSmileExcluded, nPout=nPout, nPoutExcluded=nPoutExcluded,
                      meanSmile=meanSmile, meanPout=meanPout, meanRating=meanRating, sdSmile=sdSmile, sdPout=sdPout, tValue=tValue, pValue=pValue, df=df,
                      BFplus0=BFplus0, BFr0=BFr0, d=d, dSE=dSE, dLowerCI=dLowerCI, dUpperCI=dUpperCI, g=g, gSE=gSE, gLowerCI=gLowerCI, gUpperCI=gUpperCI,
                      rawMeanDiff=rawMeanDiff, rawMeanDiffSE=rawMeanDiffSE, dataRaw=dataRaw, dataAnalysis=dCleaned4, filename=filename)
  
  return(studyResult)
  
}

extractValues <- function(results, exclusionStatus, what) {
  #
  # function that extracts values from a results list where each list element corresponds to the output of facialFeedbackAnalysis(), i.e., each list element
  # corresponds to the analysis of one study
  #
  # input:  results: results list where each list element corresponds to the output of facialFeedbackAnalysis()
  #         exclusionStatus: either "withExclusion" or "withoutExclusion", indicates which analysis to use
  #         what: the name of the list element to extract, possible are
  #
  #           "ratingsSmile": participants' cartoon ratings in the smile condition
  #           "ratingsPout": participants' cartoon ratings in the pout condition
  #           "nSmile": number of participants in the smile condition
  #           "nSmileExcluded": number of participants in the smile condition that were excluded
  #           "nPout": number of participants in the pout condition
  #           "nPoutExcluded": number of participants in the pout condition that were excluded
  #           "meanSmile": mean cartoon rating in the smile condition
  #           "meanPout": mean cartoon rating in the pout condition
  #           "meanRating": overall mean rating of participants (i.e., based on participants in both conditions)
  #           "sdSmile": standard deviation of the cartoon ratings in the smile condition
  #           "sdPout": standard deviation of the cartoon ratings in the pout condition
  #           "tValue": t-value of a one-sided independent samples t-test (equal variances assumed) which tests the hypothesis that the smile condition has a larger mean rating
  #           "pValue": p-value of the t-test
  #           "df": degrees of freedom of the t-test
  #           "BFplus0": independent samples t-test Bayes factor in favor of the hypothesis that the smile condition has a larger mean rating ( cauchy prior width = 1/sqrt(2) )
  #           "BFr0": replication Bayes factor as proposed in Verhagen & Wagenmakers (2014)
  #           "d": Cohen's d effect size
  #           "dSE": standard error of d
  #           "dLowerCI": lower bound of a 95% confidence interval for d
  #           "dUpperCI": upper bound of a 95% confidence interval for d
  #           "g": Hedges' g effect size
  #           "gSE": standard error of g
  #           "gLowerCI": lower bound of a 95% confidence interval for g
  #           "gUpperCI": upper bound of a 95% confidence interval for g
  #           "rawMeanDiff": difference between the mean of smile and pout condition
  #           "rawMeanDiffSE": standard error of difference between the mean of smile and pout condition
  #
  # output: numeric vector with the requested value for each study
  #
  
  values <- numeric(length(results))
  
  for (i in seq_along(results)) {
    
    values[i] <- results[[i]][[exclusionStatus]][[what]]
    
  }
  
  return(values)
  
}

#@@ ANALYZE ALL STUDIES @@##

# assumes that in the working directory there is a folder called "Data" which contains the .csv files with the data
# (only works if the data folder only! contains the to be analyzed .csv files and no other files!)

files <- list.files("Data")
studyIDs <- character(length(files))
results <- list()
nStudies <- length(files)

for (i in seq_along(files)) {
  
  studyIDs[i] <- paste0(i, "_", strsplit(files[i], split = ".csv")[[1]][1])
  filename <- paste0("Data/", files[i])
  results[[studyIDs[i]]][["withExclusion"]] <- facialFeedbackAnalysis(filename, excludeParticipants=TRUE)
  results[[studyIDs[i]]][["withoutExclusion"]] <- facialFeedbackAnalysis(filename, excludeParticipants=FALSE)
  
}

#@@ EXTRACT VALUES @@##

# with exclusion

meanSmileEx <- extractValues(what="meanSmile", results=results, exclusionStatus="withExclusion")
meanPoutEx <- extractValues(what="meanPout", results=results, exclusionStatus="withExclusion")
meanRatingEx <- extractValues(what="meanRating", results=results, exclusionStatus="withExclusion")
sdSmileEx <- extractValues(what="sdSmile", results=results, exclusionStatus="withExclusion")
sdPoutEx <- extractValues(what="sdPout", results=results, exclusionStatus="withExclusion")
nSmileEx <- extractValues(what="nSmile", results=results, exclusionStatus="withExclusion")
nSmileExcludedEx <- extractValues(what="nSmileExcluded", results=results, exclusionStatus="withExclusion")
nPoutEx <- extractValues(what="nPout", results=results, exclusionStatus="withExclusion")
nPoutExcludedEx <- extractValues(what="nPoutExcluded", results=results, exclusionStatus="withExclusion")

tValueEx <- extractValues(what="tValue", results=results, exclusionStatus="withExclusion")
dfEx <- extractValues(what="df", results=results, exclusionStatus="withExclusion")
pValueEx <- extractValues(what="pValue", results=results, exclusionStatus="withExclusion")

BFplus0Ex <- extractValues(what="BFplus0", results=results, exclusionStatus="withExclusion")
BFr0Ex <- extractValues(what="BFr0", results=results, exclusionStatus="withExclusion")

gEx <- extractValues(what="g", results=results, exclusionStatus="withExclusion")
gSEEx <- extractValues(what="gSE", results=results, exclusionStatus="withExclusion")
gLowerCIEx <- extractValues(what="gLowerCI", results=results, exclusionStatus="withExclusion")
gUpperCIEx <- extractValues(what="gUpperCI", results=results, exclusionStatus="withExclusion") 

dEx <- extractValues(what="d", results=results, exclusionStatus="withExclusion")
dSEEx <- extractValues(what="dSE", results=results, exclusionStatus="withExclusion")
dLowerCIEx <- extractValues(what="dLowerCI", results=results, exclusionStatus="withExclusion")
dUpperCIEx <- extractValues(what="dUpperCI", results=results, exclusionStatus="withExclusion") 

rawMeanDiffEx <- extractValues(what="rawMeanDiff", results=results, exclusionStatus="withExclusion")
rawMeanDiffSEEx <- extractValues(what="rawMeanDiffSE", results=results, exclusionStatus="withExclusion")

#@@ SAVE RESULTS TO FILE @@##

resultDataFrame <- data.frame(site=studyIDs, treat=meanSmileEx, control=meanPoutEx, sdtreat=sdSmileEx, sdcontrol=sdPoutEx, 
                              ntreat=nSmileEx, ntreatexcl=nSmileExcludedEx, ncontrol=nPoutEx, ncontrolexcl=nPoutExcludedEx, 
                              md=rawMeanDiffEx, sdmd=rawMeanDiffSEEx, t=tValueEx, df=dfEx, pvalue=pValueEx, g=gEx, SEg=gSEEx,
                              d=dEx, SEd = dSEEx)
original <- c("original", 5.14, 4.32, NA, NA, 92, 1, 92, 1, .82, .4432, 1.85, 89, .03, 0.0806, .208, .081, .2097)
resultDataFrame$site <- as.character(resultDataFrame$site)
resultDataFrame <- rbind(resultDataFrame, original)
resultDataFrame$SEg <- as.numeric(as.character(resultDataFrame$SEg))
resultDataFrame$SEd <- as.numeric(as.character(resultDataFrame$SEd))

resultDataFrame <- mutate(resultDataFrame,
                          vg = SEg^2,
                          vd = SEd^2,
                          es = "md",
                          replicated = 0,
                          experiment = "Wagenmakers")
resultDataFrame$site = gsub("[0-9]+_", "", gsub("_Data", "", resultDataFrame$site))

write.csv(resultDataFrame, "../../rrr_wagenmakers.csv", row.names = F)
