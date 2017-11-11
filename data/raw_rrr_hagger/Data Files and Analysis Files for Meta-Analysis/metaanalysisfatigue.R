AnalysisDir <- 'E:/meta'
### Set Working Directory
setwd(AnalysisDir)

#import the metafor meta-analysis R library. You will need to install it first.
require("metafor")

# Import the data
fatdat <- read.csv("fatigue.csv")
fatdat <- effdat[order(fatdat$study),]
### random-effects model meta-analysis 
# modeled after http://www.metafor-project.org/doku.php/tips:assembling_data_smd
effectSizesAll<- escalc(measure="SMD", #standardized mean difference
                     m1i= Ego.Depletion.Mean, m2i= Control.Mean,
                     sd1i=Ego.Depletion.Std.Dev, sd2i= Control.Std.Dev,
                     n1i=Ego.Depletion.Sample.size, n2i=Control.Sample.size,
                     data= fatdat)
#(fatdat$Ego.Depletion.Mean - fatdat$Control.Mean) / fatdat$Ego.Depletion.Std.Dev

res <- rma(data=effectSizesAll, yi,vi,  method="REML", slab=paste(Study.name))

#res <- rma(measure="RD", method="REML")

#For a standard random-effects model, we need to add random effects for the trials, which can be done with:
#  rma.mv(yi, vi, random = ~ 1 | trial, data=dat)

###############################################################
### FOREST PLOT WITH META-ANALYSIS ACROSS REPLICATION STUDIES
###############################################################

### Create the forest plot. 
### Rows specify the rows in which the subsets of studies appear.
### addfit=false eliminates the meta-analytic result for now. We'll first not calculate that
### so we can create a forest plot that includes all the studies, even the original, and
### then we'll calculate the meta-analysis based just on the RRR studies.

### decrease margins so the full space is used and set the font size for the forest plot
par(mar=c(4,4,1,2))
par(cex=.8, font=1)

#Create a forest plot that includes the original study
forest(res,
       ilab=cbind(fatdat$Control.Mean, fatdat$Ego.Depletion.Mean),
       addfit=FALSE, #Don't include polygon because that would include original
       ilab.xpos=c(-1.2,-1.5), #ylim=c(-1, 47) 
       #xlim=c(-10,7)
       #xlim=c(-12,5)
       #xlim=c(-12,8)
       #xlim=c(-5,2.5)
      xlim=c(-3,2.2)
)
par("usr") #reveals xlims? http://r.789695.n4.nabble.com/tweaking-forest-plot-metafor-package-td4636818.html


### add a polygon to the Forest plot showing the meta-analytic effect of the replication studies
addpoly(res, atransf=FALSE, row=0, cex=1.0, mlab="Meta-analytic effect")

abline(h=0.45)

