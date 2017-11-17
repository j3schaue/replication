###-----------------------------------------------------------###
###-----------------------------------------------------------###
### Table 1S from Camerer
###  Reports effect sizes (r) and sample sizes (n) per lab
###  Cobbled from PDF in paper.
###-----------------------------------------------------------###
###-----------------------------------------------------------###

# set the working directory to src so we can use relative paths
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../raw/rpe")

source("S1pvalsetc.R") # lists of p-values from the paper
orig = read.csv("S1origes.csv") # partial table copied from paper PDF
reps = read.csv("S1repes.csv") # partial table copied from paper PDF
exps = read.csv("S1experiment.csv") # partial table copied from paper PDF
repdet = read.csv("S1replicationdet.csv") # partial table copied from paper PDF

#---------------------------------------------#
## Combine information to re-create table S1
#---------------------------------------------#
s1 = do.call(cbind, list(exps, orig, reps, repdet))

# save S1
write.csv(s1, "../camerer_S1.csv", row.names=F)

#---------------------------------------------#
## Modify table S1 for analysis
#---------------------------------------------#
library(tidyr); library(dplyr)

#---Melt S1 so that each experiment has its own row
sizes = s1 %>% select(-reles, -r, -rrep) %>% gather(replicate, n, c(n,nrep))
sizes$replicate = as.integer(sizes$replicate == "nrep")
corrs = s1 %>% select(-replicated, -reles, -n, -nrep) %>% gather(replicate, r, c(r,rrep))
corrs$replicate = as.integer(corrs$replicate == "rrep")
df = left_join(sizes, corrs)  # melted df
df$n = as.numeric(as.character(df$n))

#---Transform correlations 
df$z = 0.5*log((1 + df$r)/(1 - df$r)) # Fisher transfrom
df$vz = 1/(df$n - 3)
head(df)
df %>% group_by(experiment) %>% summarize(mn = mean(replicated=='No'))

#---Code replicate names
grabName = function(string){
  # Function that takes names of experiments in THIS df and returns just the 
  # first word (i.e., first author's last name)
  # This only works since these are unique for each experiment!
  nm = strsplit(string, " ")[[1]][1]
  if(nm == "de"){ # conditional for de Clippel
    nm = paste0(strsplit(string, " ")[[1]][1], strsplit(string, " ")[[1]][2])
  }
  return(nm)
}

nms = sapply(as.character(df$experiment), grabName) # list of names
reps = sapply(df$replicate, FUN=function(x) ifelse(x==0, "_orig", "_rep")) # indicate replicates
df$site = paste0(nms, reps) # site = name_replicate format

# rename experiment
df$experiment = gsub("[[:blank:]]*\\([[:alnum:]]*[[:blank:]]+[[:alnum:]]*\\)", "", df$experiment)

# denote ES type
df$es = "z"

# standardize replcated column
df$replicated = as.integer(df$replicated == 'Yes')

str(df)
df

write.csv(df, "../../rpe.csv", row.names=F)

# write.csv(df %>% select(experiment, site, es, z, vz, replicated), "../camerer.csv", row.names=F)

