###-----------------------------------------------------------###
###-----------------------------------------------------------###
### Table 1S from Camerer
###  Reports effect sizes (r) and sample sizes (n) per lab
###  Cobbled from PDF in paper.
###-----------------------------------------------------------###
###-----------------------------------------------------------###

# setwd("~/Documents/replication/replication/data/raw_camerer/") #path on Jake's machine
setwd("./data/raw_camerer") # relative path
list.files()


source("S1pvalsetc.R") # lists of p-values from the paper
orig = read.csv("S1origes.csv") # partial table copied from paper PDF
reps = read.csv("S1repes.csv") # partial table copied from paper PDF
exps = read.csv("S1experiment.csv") # partial table copied from paper PDF
repdet = read.csv("S1replicationdet.csv") # partial table copied from paper PDF

# Combine information to re-create table S1
s1 = do.call(cbind, list(exps, orig, reps, repdet))

# save S1
write.csv(s1, "../camerer_S1.csv", row.names=F)


# Modify table for analysis
library(tidyr); library(dplyr)

# melt S1 so that each experiment has its own row
sizes = s1 %>% select(-replicated, -reles, -r, -rrep) %>% gather(replicate, n, c(n,nrep))
sizes$replicate = as.integer(sizes$replicate == "nrep")
corrs = s1 %>% select(-replicated, -reles, -n, -nrep) %>% gather(replicate, r, c(r,rrep))
corrs$replicate = as.integer(corrs$replicate == "rrep")
df = left_join(sizes, corrs)  # melted df

# transform correlations 
df$z = 0.5*log((1 + df$r)/(1 - df$r))
df$vz = 1/(df$n - 3)
head(df)

# Code replicate names
grabName = function(string){
  nm = strsplit(string, " ")[[1]][1]
  if(nm == "de"){
    nm = paste0(strsplit(string, " ")[[1]][1], strsplit(string, " ")[[1]][2])
  }
  return(nm)
}
nms = sapply(as.character(df$experiment), grabName)
reps = sapply(df$replicate, FUN=function(x) ifelse(x==0, "_orig", "_rep"))
df$site = paste0(nms, reps)

# rename experiment
df$experiment = gsub("[[:blank:]]*\\([[:alnum:]]*[[:blank:]]+[[:alnum:]]*\\)", "", df$experiment)

# denote ES type
df$es = "r"

str(df)

write.csv(df, "camerer_full.csv")

write.csv(df %>% select(experiment, site, es, z, vz), "../camerer.csv", row.names=F)
