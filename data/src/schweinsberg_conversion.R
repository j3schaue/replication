###-----------------------------------------------------------###
###-----------------------------------------------------------###
### PPIR files
###-----------------------------------------------------------###
###-----------------------------------------------------------###
library(foreign); library(dplyr); library(ggplot2)
tmp = read.spss("1 effectsizes.graphdata 2015 08 31 1431.sav")

names(tmp)
sapply(tmp, length)
df = data.frame(tmp)
df$Study = trimws(df$Study)
df$Sample = trimws(df$Sample)
write.csv(df, "Shweinsberg_fig1.csv", row.names=F)
unique(df$Study)
head(df)
ggplot(df) + aes(x = EffectSize, y=Study) + geom_point() #+ 
  geom_point(aes(x=OriginalES, y=Study), color='red')
