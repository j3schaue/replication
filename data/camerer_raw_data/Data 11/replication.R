########## 
# Friedman and Oprea (2012) Repplication
# 
# Colin Camerer, Taisuke Imai, Gideon Nave
# Last update: October 12, 2015
########## 

setwd("SET WORKING DIRECTORY HERE")

# Import data
D           <- read.table("profiles_rep.csv",header=FALSE,sep=",")
colnames(D) <- c("session", "inst", "Period", "uniquesub1", "uniquesub2", "phase", "Tick", "switch", "time.share", "P1S", "P2S", "y", "x", "subp", "group", "duplicate", "treat")

# Treatment information
D$inst <- ifelse(D$inst==1, "C", "G")

# Condition information
D$treat <- ifelse(D$treat==1,"Easy", ifelse(D$treat==2, "Mix-a", ifelse(D$treat==3, "Mix-b", "Hard")))

# Strategy profile
D$profile <- ifelse(D$P1S==1 & D$P2S==1, "coop", ifelse(D$P1S==0 & D$P2S==0, "def", "off"))

# Look at periods 13-32 (behavior stabilized)
per <- 12
D   <- D[D$Period > per,]

# Calculate time weighted rates of cooperation, defection and ST (called "off")
D$time.share <- D$time.share/60
D$Tick       <- D$Tick/D$subp

D$coop0 <- (D$P1S==1 & D$P2S==1)*D$time.share
D$def0  <- (D$P1S==0 & D$P2S==0)*D$time.share
D$off0  <- ((D$P1S==1 & D$P2S==0) | (D$P1S==0 & D$P2S==1))*D$time.share

max(D$coop+D$def+D$off)

D$coop <- ave(D$coop0,D$group,FUN=sum)
D$def  <- ave(D$def0,D$group,FUN=sum)
D$off  <- ave(D$off0,D$group,FUN=sum)

# Calculate final times of cooperation
D$coop1               <- (D$profile=="coop")*1
D$coop2               <- D$coop1*(D$Tick+D$time.share)
D[D$inst=="G",]$coop2 <- D[D$inst=="G",]$coop1*D[D$inst=="G",]$Tick
D$lastcoop            <- ave(D$coop2,D$group,FUN=max)
D$cut                 <- D$lastcoop

D$mintick <- ave(D$Tick,D$group,FUN=min)
dat       <- D[D$Tick==D$mintick,]

dat$xhigh <- (dat$x==18)*1
dat$yhigh <- (dat$y==8)*1

# Rename 
panel0 <- aggregate.data.frame(list(coop=dat$coop,cut=dat$cut),list(id=dat$uniquesub1,inst=dat$inst),FUN=median)
panel  <- aggregate.data.frame(list(coop=dat$coop,cut=dat$cut),list(id=dat$uniquesub1,inst=dat$inst,treat=dat$treat,x=dat$x,y=dat$y),FUN=median)

panel$xhigh <- (panel$x==18)*1
panel$yhigh <- (panel$y==8)*1

tapply(panel0$cut,panel0$inst,FUN=median)
tapply(panel$cut,list(inst=panel$inst,treat=panel$treat),FUN=median)

# Test Hypothesis of Difference Across Treatments
wilcox.test(panel0[panel0$inst=="C",]$coop,panel0[panel0$inst=="G",]$coop)

