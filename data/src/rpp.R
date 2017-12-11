###-----------------------------------------------------------###
###-----------------------------------------------------------###
### RPP files
###-----------------------------------------------------------###
###-----------------------------------------------------------###

library(dplyr); library(tidyr)

# set the working directory to src so we can use relative paths
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../raw/")

rp = read.csv("./rpp.csv")
unique(rp$Study.Num)

#---Select the relevant data
arp = rp %>% 
  select(authors=Authors..O., num=Study.Num, experiment=Study.Title..O., 
         r=T_r..O., n=N..O., rrep=T_r..R., nrep=N..R., 
         pvalo = T_pval_USE..O., pvalr = T_pval_USE..R., 
         stat=T_Test.Statistic..R., df1=T_df1..O.,
         replicated=Replicate..R., cirep=O.within.CI.R, 
         meta=Meta.analysis.significant) %>%
  mutate(replicated=tolower(as.character(replicated))) %>%
  filter((stat=='F' & df1==1) | stat=='t' | stat=='r') # get the meta-analytic subset

dim(arp)

#---Melt dataframe so that each experiment has its own row
# df of experiments and sample sizes
sizes = arp %>% select(-r, -rrep) %>% gather(replicate, n, c(n,nrep))
sizes$replicate = as.integer(sizes$replicate == "nrep") # denote replicates

# df of experiments and effect sizes
corrs = arp %>% select(-n, -nrep) %>% gather(replicate, r, c(r,rrep))
corrs$replicate = as.integer(corrs$replicate == "rrep") # denote replicates

# merge sample sizes and effect sizes
df = left_join(sizes, corrs)  # melted df
df$n = as.numeric(as.character(df$n)) # convert sample sizes to numeric

#---Transform correlations 
df$z = 0.5*log((1 + df$r)/(1 - df$r)) # Fisher transform
df$vz = 1/(df$n - 3)
head(df)

#---Clean up author names
authors = lapply(df$authors, FUN=function(x) 
  strsplit(as.character(gsub('\xd5', 'i', x)), ",")[[1]][1])

which(duplicated(authors))
dups = which(authors %in% authors[c(63, 65)])
df[dups, c('authors', 'num')]
authors[c(17, 90)]  = paste(authors[17], 1)
authors[c(63, 136)]  = paste(authors[63], 2)
authors[c(64, 137)]  = paste(authors[64], 1)
authors[c(65, 138)]  = paste(authors[65], 2)

df$authors = unlist(authors)

#---Mark experiment name
out = df %>%
  rename(exp_name = experiment) %>% 
  rename(experiment = authors) %>% 
  mutate(replicate = ifelse(replicate==1, "_rep", "orig")) %>%
  unite(site, experiment, replicate, sep="_") %>% 
  cbind(., data.frame(experiment=df$authors)) %>%
  select(experiment, site, n, r, z, vz, replicated, cirep, 
         meta, stat, pvalo, pvalr, exp_name) %>%
  mutate(replicated = as.integer(replicated == 'yes'))
out$site = gsub("__", "_", out$site)

head(out)
dim(out)

#---denote ES type
out$es = "z"

#---Write to file
write.csv(out, "../rpp.csv")

