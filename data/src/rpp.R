###-----------------------------------------------------------###
###-----------------------------------------------------------###
### RPP files
###-----------------------------------------------------------###
###-----------------------------------------------------------###

library(dplyr); library(tidyr)

setwd("~/Documents/replication/replication/data/raw") #jakes path
# setwd("./data") # relative path

rp = read.csv("./rpp.csv")
unique(rp$Study.Num)

#---Select the relevant data
arp = rp %>% 
  select(num=Study.Num, experiment=Study.Title..O., 
         r=T_r..O., n=N..O., rrep=T_r..R., nrep=N..R., pvalo = T_pval_USE..O., pvalr = T_pval_USE..R., stat=T_Test.Statistic..R., df1=T_df1..O.,
         rep1=Replicate..R., cirep=O.within.CI.R, meta=Meta.analysis.significant) %>%
  filter((stat=='F' & df1==1) | stat=='t' | stat=='r') # get the meta-analytic subset

dim(arp)

#---Melt dataframe so that each experiment has its own row
sizes = arp %>% select(-r, -rrep) %>% gather(replicate, n, c(n,nrep))
sizes$replicate = as.integer(sizes$replicate == "nrep") # denote replicates
corrs = arp %>% select(-n, -nrep) %>% gather(replicate, r, c(r,rrep))
corrs$replicate = as.integer(corrs$replicate == "rrep") # denote replicates
df = left_join(sizes, corrs)  # melted df
df$n = as.numeric(as.character(df$n)) # convert sample sizes to numeric

#---Transform correlations 
df$z = 0.5*log((1 + df$r)/(1 - df$r)) # Fisher transform
df$vz = 1/(df$n - 3)
head(df)

#---Mark experiment name
out = df %>% rename(exp_name = experiment) %>% rename(experiment = num) %>% 
  unite(site, experiment, replicate, sep="_") %>% cbind(experiment=df$num) 

head(out)
dim(out)

#---denote ES type
out$es = "r"

write.csv(out[, c('experiment', 'exp_name', 'site', 'es', 'stat', 'df1', 'n', 'r', 'z', 'vz', 'rep1', 'cirep', 'meta', 'pvalo', 'pvalr')], 
          "../rpp.csv")

