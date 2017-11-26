## N.B.: The scales used for producing the figures may change
## for the analysis of the real data, however, the core
## analysis will not change. 

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

##@@ load Josine Verhagen's Replication Bayes Factor Functions @@##

source("Repfunctionspack.R")

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

# without exclusion

meanSmileAll <- extractValues(what="meanSmile", results=results, exclusionStatus="withoutExclusion")
meanPoutAll <- extractValues(what="meanPout", results=results, exclusionStatus="withoutExclusion")
meanRatingAll <- extractValues(what="meanRating", results=results, exclusionStatus="withoutExclusion")
sdSmileAll <- extractValues(what="sdSmile", results=results, exclusionStatus="withoutExclusion")
sdPoutAll <- extractValues(what="sdPout", results=results, exclusionStatus="withoutExclusion")
nSmileAll <- extractValues(what="nSmile", results=results, exclusionStatus="withoutExclusion")
nSmileExcludedAll <- extractValues(what="nSmileExcluded", results=results, exclusionStatus="withoutExclusion")
nPoutAll <- extractValues(what="nPout", results=results, exclusionStatus="withoutExclusion")
nPoutExcludedAll <- extractValues(what="nPoutExcluded", results=results, exclusionStatus="withoutExclusion")

tValueAll <- extractValues(what="tValue", results=results, exclusionStatus="withoutExclusion")
dfAll <- extractValues(what="df", results=results, exclusionStatus="withoutExclusion")
pValueAll <- extractValues(what="pValue", results=results, exclusionStatus="withoutExclusion")

BFplus0All <- extractValues(what="BFplus0", results=results, exclusionStatus="withoutExclusion")
BFr0All <- extractValues(what="BFr0", results=results, exclusionStatus="withoutExclusion")

gAll <- extractValues(what="g", results=results, exclusionStatus="withoutExclusion")
gSEAll <- extractValues(what="gSE", results=results, exclusionStatus="withoutExclusion")
gLowerCIAll <- extractValues(what="gLowerCI", results=results, exclusionStatus="withoutExclusion")
gUpperCIAll <- extractValues(what="gUpperCI", results=results, exclusionStatus="withoutExclusion") 

dAll <- extractValues(what="d", results=results, exclusionStatus="withoutExclusion")
dSEAll <- extractValues(what="dSE", results=results, exclusionStatus="withoutExclusion")
dLowerCIAll <- extractValues(what="dLowerCI", results=results, exclusionStatus="withoutExclusion")
dUpperCIAll <- extractValues(what="dUpperCI", results=results, exclusionStatus="withoutExclusion") 

rawMeanDiffAll <- extractValues(what="rawMeanDiff", results=results, exclusionStatus="withoutExclusion")
rawMeanDiffSEAll <- extractValues(what="rawMeanDiffSE", results=results, exclusionStatus="withoutExclusion")

#@@ SAVE RESULTS TO FILE @@##

resultDataFrame <- data.frame(studyIDs=studyIDs, meanSmileEx=round(meanSmileEx, 3), meanPoutEx=round(meanPoutEx, 3), meanRatingEx=round(meanRatingEx, 3),
                              sdSmileEx=round(sdSmileEx, 3), sdPoutEx=round(sdPoutEx, 3), nSmileEx=nSmileEx, nSmileExcluded=nSmileExcludedEx, nPoutEx=nPoutEx,
                              nPoutExcluded=nPoutExcludedEx, tValueEx=round(tValueEx, 3), dfEx=round(dfEx, 3), pValueEx=round(pValueEx, 3),
                              BFplus0Ex=round(BFplus0Ex, 3), BFr0Ex=round(BFr0Ex, 3), gEx=round(gEx, 3), gSEEx=round(gSEEx, 3), gLowerCIEx=round(gLowerCIEx, 3),
                              gUpperCIEx=round(gUpperCIEx, 3), dEx=round(dEx, 3), dLowerCIEx=round(dLowerCIEx, 3), dUpperCIEx=round(dUpperCIEx, 3),
                              rawMeanDiffEx=round(rawMeanDiffEx, 3), rawMeanDiffSEEx=round(rawMeanDiffSEEx, 3), meanSmileAll=round(meanSmileAll, 3), meanPoutAll=round(meanPoutAll, 3),
                              meanRatingAll=round(meanRatingAll, 3), sdSmileAll=round(sdSmileAll, 3), sdPoutAll=round(sdPoutAll, 3), nSmileAll=nSmileAll, nPoutAll=nPoutAll,
                              tValueAll=round(tValueAll, 3), dfAll=round(dfAll, 3), pValueAll=round(pValueAll, 3), BFplus0All=round(BFplus0All, 3),
                              BFr0All=round(BFr0All, 3), gAll=round(gAll, 3), gSEAll=round(gSEAll, 3), gLowerCIAll=round(gLowerCIAll, 3), gUpperCIAll=round(gUpperCIAll, 3), dAll=round(dAll, 3),
                              dLowerCIAll=round(dLowerCIAll, 3), dUpperCIAll=round(dUpperCIAll, 3), rawMeanDiffAll=round(rawMeanDiffAll, 3), rawMeanDiffSEAll=round(rawMeanDiffSEAll, 3))

write.csv(resultDataFrame, file="resultsFacialFeedbackReplication.csv", row.names=FALSE, quote=FALSE)


#@@ CALCULATE META-ANALYTIC EFFECT SIZE: RAW EFFECT SIZE @@##

esEx <- rawMeanDiffEx
esSEEx <- rawMeanDiffSEEx
xlabel <- "Smile-Pout"

# random effects model
metaEx <- rma(yi = esEx, sei = esSEEx)

#@@ FOREST PLOT: RAW EFFECT SIZE @@##

# raw effect size original study
SMSes <- 5.14 - 4.32

# standard error of raw effect size original study
SMSesSE <- SMSes/1.85

### with exclusion

# trim means to two decimal digits
cMeanSmileEx <- format(meanSmileEx, digits=3)
cMeanPoutEx <- format(meanPoutEx, digits=3)

cairo_pdf("forestPlot_withExclusion_raw.pdf", width=830/72, height=530/72)

forest(x = c(SMSes, esEx), sei = c(SMSesSE, esSEEx), xlab=xlabel, cex.lab=1.4,
       ilab=cbind(c("5.14", cMeanSmileEx), c("4.32", cMeanPoutEx)),
       ilab.xpos=c(grconvertX(.18, from = "ndc", "user"),
                   grconvertX(.28, from = "ndc", "user")), cex.axis=1.1, lwd=1.4,
       rows=c(nStudies+7, (nStudies+2):3), ylim=c(-2, nStudies+11),
       slab = c("SMS Study 1\n(Original Study)", paste0("Study ", seq_len(nStudies))))

abline(h=nStudies+5, lwd=1.4)
text(grconvertX(.019, from = "ndc", "user"), nStudies+3.75, "RRR Studies", cex=1.2, pos = 4)
text(grconvertX(.053, from = "ndc", "user"), nStudies+10, "Study", cex=1.2)
text(grconvertX(.18, from = "ndc", "user"), nStudies+10, "Smile", cex=1.2)
text(grconvertX(.28, from = "ndc", "user"), nStudies+10, "Pout", cex=1.2)
text(grconvertX(.9, from = "ndc", "user"), nStudies+10, paste0(xlabel, " [95% CI]"), cex=1.2)

abline(h=1, lwd=1.4)
addpoly(metaEx, atransf=FALSE, row=-1, cex=1.3, mlab="Meta-Analytic Effect Size:")

dev.off()


#@@ CALCULATE META-ANALYTIC EFFECT SIZE: STANDARDIZED EFFECT SIZE @@##

esEx <- dEx
esSEEx <- dSEEx
xlabel <- "Cohen's d"

# standardized effect size original study
silent <- capture.output( effectSize <- tes(t=1.85, n.1=32, n.2=32) ) # capture.output allows one to suppress printing to console
SMSd <- effectSize$d

# standard error of standardized effect size original study
SMSdSE <- sqrt(effectSize$var.d)

# random effects model
metaEx <- rma(yi = esEx, sei = esSEEx)

#@@ FOREST PLOT: STANDARDIZED EFFECT SIZE @@##

### with exclusion

# trim means to two decimal digits
cMeanSmileEx <- format(meanSmileEx, digits=3)
cMeanPoutEx <- format(meanPoutEx, digits=3)

cairo_pdf("forestPlot_withExclusion_standardized.pdf", width=830/72, height=530/72)

forest(x = c(SMSd, esEx), sei = c(SMSdSE, esSEEx), xlab=xlabel, cex.lab=1.4,
       ilab=cbind(c("5.14", cMeanSmileEx), c("4.32", cMeanPoutEx)),
       ilab.xpos=c(grconvertX(.18, from = "ndc", "user"),
                   grconvertX(.28, from = "ndc", "user")), cex.axis=1.1, lwd=1.4,
       rows=c(nStudies+7, (nStudies+2):3), ylim=c(-2, nStudies+11),
       slab = c("SMS Study 1\n(Original Study)", paste0("Study ", seq_len(nStudies))))

abline(h=nStudies+5, lwd=1.4)
text(grconvertX(.019, from = "ndc", "user"), nStudies+3.75, "RRR Studies", cex=1.2, pos = 4)
text(grconvertX(.053, from = "ndc", "user"), nStudies+10, "Study", cex=1.2)
text(grconvertX(.18, from = "ndc", "user"), nStudies+10, "Smile", cex=1.2)
text(grconvertX(.28, from = "ndc", "user"), nStudies+10, "Pout", cex=1.2)
text(grconvertX(.9, from = "ndc", "user"), nStudies+10, paste0(xlabel, " [95% CI]"), cex=1.2)

abline(h=1, lwd=1.4)
addpoly(metaEx, atransf=FALSE, row=-1, cex=1.3, mlab="Meta-Analytic Effect Size:")

dev.off()


#@@ 2-PANEL PLOT @@##

# packages("devtools")
# install_github("ndphillips/yarrr")
library(yarrr)
source("pirateplot_custom.R")
r <- resultDataFrame

scatterplot <- function(xVar, yVar, cexPoints = 1.3, cexXAxis = 1.3,
                        cexYAxis = 1.3, lwd = 2) {
  
  op <- par(mar = c(5, 5, 5, 2))
  low_tmp <- min(c(xVar, yVar))
  high_tmp <- max(c(xVar, yVar))
  ticks_tmp <- pretty(c(low_tmp, high_tmp))
  lim <- range(ticks_tmp)
  ticks <- pretty(lim)
  
  plot(xVar, yVar, pch = 21, bg = "grey", ylab = "", xlab = "", axes = FALSE,
       ylim = lim, xlim = lim, cex = cexPoints, asp= 1, xaxs = "i", yaxs = "i")
  lines(lim, lim, lwd = lwd)
  
  axis(1, line = .4, at = ticks, cex.axis = cexXAxis)
  axis(2, line = -.8, at = ticks, cex.axis = cexYAxis, las = 1)
  
}


cairo_pdf("2panelPlot.pdf", width = 2*400/72, height = 400/72)
op <- par(mfrow = c(1, 2))

### Panel 1: Pirate Plot
op1 <- par(mar = c(5, 6, 5, 2))
dummy <- factor(rep("", 20))
data <- data.frame(rawMeanDiffs = r$rawMeanDiffEx, dummy = dummy)
ylim <- c(-.5, 2)
p <- pirateplot_custom(rawMeanDiffs ~ dummy, data,
                       xlab = "",
                       ylab = "",
                       ylim = ylim,
                       bean.o = .8, point.o = .8,
                       pal = "darkblue",
                       yaxt = "n",
                       bar.o = 0,
                       bty = "n")
lines(p, rep(0, 2), lwd = 2, lty = 2)
mtext("Mean Difference (Smile - Pout)", 2, line = 3.5, cex = 1.3)
axis(2, las = 1, line = 0.2, cex.axis = 1.3)
par(op1)

### Panel 2: Scatterplot
means_smile <- r$meanSmileEx
means_pout <- r$meanPoutEx
scatterplot(means_pout, means_smile)
mtext("Mean in Pout Condition", 1, line = 3, cex = 1.3)
mtext("Mean in Smile Condition", 2, line = 2.8, cex = 1.3)

par(op)
dev.off()


#@@ Rating Difference vs. Mean Rating Plot @@##

scatterplot2 <- function(xVar, yVar, cexPoints = 1.3, cexXAxis = 1.2,
                         lwd = 2, specifiedXlim = NULL,
                         specifiedYlim = NULL, cexYAxis = 1.2) {
  
  op <- par(mar = c(5, 5, 2, 2))
  
  if ( ! is.null(specifiedXlim)) {
    xticks <- pretty(specifiedXlim)
  } else {
    xticks <- pretty(xVar)
  }
  
  if ( ! is.null(specifiedYlim)) {
    yticks <- pretty(specifiedYlim)
  } else {
    yticks <- pretty(yVar)
  }
  
  xlim <- range(xticks)
  ylim <- range(yticks)
  
  plot(xVar, yVar, pch = 21, bg = "grey", ylab = "", xlab = "", axes = FALSE,
       ylim = ylim, xlim = xlim, cex = cexPoints)
  axis(2, line = -.8, at = yticks, cex.axis = cexYAxis, las = 1)
  
}

cexXAxis <- 1.2

cairo_pdf("DifferenceVsMeanRatingPlot.pdf", width = 510/72, height = 400/72)
scatterplot2(r$meanRatingEx, r$rawMeanDiffEx,
             specifiedXlim = c(0, 9))
axis(1, line = .4, at = 0:9, cex.axis = cexXAxis)
mtext("Average Rating", 1, cex = 1.55, line = 2.7)
mtext("Mean Difference (Smile - Pout)", 2, cex = 1.55, line = 2.3)
dev.off()
