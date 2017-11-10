###-----------------------------------------------------------###
###-----------------------------------------------------------###
### RPP files
###-----------------------------------------------------------###
###-----------------------------------------------------------###

library(dplyr); library(tidyr)

# setwd("~/Documents/replication/replication/data/") #jakes path
setwd("./data") # relative path

rp = read.csv("raw_rpp.csv")
unique(rp$Study.Num)

# relevant data
arp = rp %>% 
  select(num=Study.Num, experiment=Study.Title..O., 
         r=T_r..O., n=N..O., rrep=T_r..R., nrep=N..R.) %>%
  filter(!is.na(r) & !is.na(n) & !is.na(rrep) & !is.na(nrep))

# melt so that each experiment has its own row
sizes = arp %>% select(-r, -rrep) %>% gather(replicate, n, c(n,nrep))
sizes$replicate = as.integer(sizes$replicate == "nrep")
corrs = arp %>% select(-n, -nrep) %>% gather(replicate, r, c(r,rrep))
corrs$replicate = as.integer(corrs$replicate == "rrep")
df = left_join(sizes, corrs)  # melted df
df$n = as.numeric(as.character(df$n))

# transform correlations 
df$z = 0.5*log((1 + df$r)/(1 - df$r))
df$vz = 1/(df$n - 3)
head(df)

# experiment name
out = df %>% rename(exp_name = experiment) %>% rename(experiment = num) %>% 
  unite(site, experiment, replicate, sep="_") %>% cbind(experiment=df$num) 

write.csv(out, "rpp_full.csv")

# denote ES type
out$es = "r"

out = out[c("experiment", "site", "es", "z", "vz")]

write.csv(out, "rpp.csv")
