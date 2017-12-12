
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(dplyr)


#experiment, n, r, replicated
ie <- c("Intuitive Economics", 225, .39, 1)
ie <- rbind.data.frame(ie)
names(ie) <- c("experiment", "n", "r", "replicated")
ie$n <- as.numeric(as.character(ie$n))
ie$r <- as.numeric(as.character(ie$r))

ie <- mutate(ie,
             site = "original",
             es = "smd",
             vr = ((1-r^2)^2)/(n-1), 
             d = 2*r/sqrt(1-r^2), 
             vd = 4*vr/((1-r^2)^3),
             g = d*(1 - 3/(4*(n - 2) - 1)),
             vg = vd*((1 - 3/(4*(n - 2) - 1))^2))
ie <- select(ie, experiment, site, n, d, vd, g, vg, es, replicated)

#Experiment, n, n1, n2, es1, es2, sd1, sd2, replicated
pog <- c("Presumption of Guilt", 79, 39, 40, 3.97, 3.93, 1.27, 1.42, 0)
mi <- c("Moral Inversion", 58, 29, 29, 4.33, 3.31, .9, 1.54, 1)
mc <- c("Moral Cliff", 107, 53, 54, 5.07, 4.14, 1.36, 1.26, 1)
hscompany <- c("Higher Standards - Company", 90, 45, 45, 3.94, 3.47, 1.25, 1.47, 1)
hscharity <- c("Higher Standards - Charity", 75, 38, 37, 4.25, 3.03, 1.29, 1.36, 0)
chp <- c("Cold-Hearted Prosociality", 80, 40, 40, 4.56, 2.04, .93, 1.27, 1) #weird analysis
bih <- c("Burn in Hell", 154, 77, 77, .42, .34, .3, .29, 1)
bt <- c("Bad Tipper", 77, 38, 39, 4.41, 3.57, 1.27, 1.35, 1)
bai <- c("Belief-Act Inconsistency", 126, 63, 63, -.92, -1.58, 1.72, 1.81, 1)

df <- rbind.data.frame(pog,mi,mc,hscompany,hscharity,chp,bih,bt,bai)
names(df) <- c("experiment", "n", "n1", "n2", "es1", "es2", "sd1", "sd2", "replicated")
df$n <- as.numeric(as.character(df$n))
df$n1 <- as.numeric(as.character(df$n1))
df$n2 <- as.numeric(as.character(df$n2))
df$es1 <- as.numeric(as.character(df$es1))
df$es2 <- as.numeric(as.character(df$es2))
df$sd1 <- as.numeric(as.character(df$sd1))
df$sd2 <- as.numeric(as.character(df$sd2))


df <- mutate(df,
             site = "original",
             es = "smd",
             d = (es1-es2)/sqrt(((n1-1)*sd1^2+(n2-1)*sd2^2)/(n1+n2-2)),
             vd = (n1+n2)/(n1*n2) + (d^2)/(2*(n1+n2)),
             J = 1 - 3/(4*(n1+n2-2) - 1),
             g = d * J,
             vg = vd * J^2)

df <- select(df, experiment, site, n, d, vd, g, vg, es, replicated)

df <- rbind(df, ie)
summary(df$vd)


#---Bigot-Misathrope
baorig = data.frame(experiment='Bigot-Misanthrope', site='original', n=46,
                    d=6.07/sqrt(45), vd=1/45, g=NA, vg=NA, es='d', replicated=1)

df = rbind(df, baorig)
str(df)
df$replicated = as.integer(as.character(df$replicated))

write.csv(df, "../raw/ppir/ppir_original.csv", row.names = FALSE)
