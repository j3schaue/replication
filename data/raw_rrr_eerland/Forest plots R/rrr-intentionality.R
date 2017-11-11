library(metafor)
data <- read.csv(file = "Intentionality.csv", header = TRUE, sep = ",", quote = "\"", dec = ".", fill = TRUE)
Authors <- data$Author
dat <- escalc(measure="MD", m1i=mimp, sd1i=sdimp, n1i=nimp, m2i=mperf, sd2i=sdperf, n2i=nperf, data=data)
res <- rma(yi, vi, data=dat)
par(mar=c(4,4,1,2))
par(cex=.7, font=1)
forest(res, alim=c(-12,6), at=c(-2, -1, 0, 1, 2), addfit=FALSE, atransf=FALSE, ilab=cbind(data$mimp,data$mperf), ilab.xpos=c(-5,-2.80), ylim=c(-1, 19), rows=c(16,14,12:2), xlab="Intentionality", mlab="Random Effects Model", slab=paste(Authors))
text(-5,18, "Imperfective", cex=.9)
text(-2.80,18.05, "Perfective", cex=.8)
text(5, 18, "Difference [95% CI]", cex=.8)
text(-11.50,18, "Study", cex=.8)
abline(h=13)
abline(h=15)
abline(h=1)
repdata <- data[3:13,]
rownames(repdata) <- NULL
dat1 <- escalc(measure="MD", m1i=mimp, sd1i=sdimp, n1i=nimp, m2i=mperf, sd2i=sdperf, n2i=nperf, data=repdata)
repres <- rma(yi, vi, data=dat1)
repres
addpoly(repres, atransf=FALSE, row=0, cex=.8, mlab="Meta-analytic effect for laboratory replications only")
