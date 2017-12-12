###-----------------------------------------------------------###
###-----------------------------------------------------------###
### PPIR files
###-----------------------------------------------------------###
###-----------------------------------------------------------###

library(dplyr); library(tidyr)

# set the working directory to src so we can use relative paths
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../raw/ppir")

#---------------------------------------------#
## Sample Sizes by University
#---------------------------------------------#
sz = read.csv("schweinsberg_samplesizes_new.csv")
names(sz)[1] = "Study"
sz$Study = as.character(sz$Study)
sz$University = as.character(sz$University)

# Clean some of the character strings to match other dfs
sz$University = gsub(" in Washington DC", "",
  gsub(", (Canada|Germany|the Netherlands|China)", "", sz$University))

#---------------------------------------------#
## University vs. individual researcher names
#---------------------------------------------#
nm = read.csv("name-uni.csv")
nm$name = as.character(nm$name)
nm$uni = as.character(nm$uni)

# standardize names and strip whitespace
nm$name = gsub("jenniferjordan", "jenjordan",
            gsub("victoriabrescoll", "toribrescoll",
              gsub("xiaominsun", "sunnysun",
                gsub("danmolden[[:alnum:]]+", "danmolden",
                  tolower(gsub("[[:blank:]]", "", nm$name))))))
# list of unique names to check against ES data file
nmsun = sort(unique(nm$name))

# Clear up inconsitencies with how names and universities were entered
nm[nm$name=="danmolden",]$uni = "Northwestern University"
nm = unique(nm) # drop replicates
nm[nm$name=="inseadsorbonnelab",] = c("french", "INSEAD, France") # inconsitencies with French institutes
nm = rbind(nm, data.frame(name=c("mturk", "ucibusiness", "ucipsychstudents"), # other inconsitencies
           uni=c("Mechanical Turk sample", "University of California Irvine", "University of California Irvine"))
)


#---------------------------------------------#
## Effect size data
#---------------------------------------------#
dat = read.csv("Shweinsberg_fig1.csv")
dat$Study = as.character(dat$Study)

# Fix cases and hyphens
dat$Study[dat$Study=="Intuitive economics"] = "Intuitive Economics"
dat$Study[dat$Study=="Cold Hearted Prosociality"] = "Cold-Hearted Prosociality"

# standardize names
dat$name = 
  gsub("annlaure", "annelauresellier",
   gsub("toribrescolli", "toribrescoll", 
    gsub("additionalfrench", "french",
     gsub("data|translation", "",
       tolower(gsub("[[:blank:]]", "", dat$Sample))))))

# check 'french' names
dat %>%
  filter(name %in% c("french", "frenchtranslation", "additionalfrenchdata"))
sz[grepl("INSEAD", sz$University), ]

#---------------------------------------------#
## check differences between DFs
#---------------------------------------------#
# names
setdiff(unique(nm$name), unique(dat$name))
setdiff(unique(dat$name), unique(nm$name))
# universities
uni1 = sort(unique(sz$University))
uni2 = sort(unique(nm$uni))
setdiff(uni1, uni2)
setdiff(uni2, uni1)

#-------------------------------------------------------------------#
## Merge dat and sz, but we need to connect the dat to universities 
#-------------------------------------------------------------------#
names(nm)[2] = "University"
tmp = left_join(dat, nm) 

# check join
tmp %>% group_by(Study, University) %>% tally() %>% filter(n>1)
tmp %>% filter(University == "University of Washington (Foster)")
apply(tmp, 2, FUN=function(x) mean(is.na(x)))

#---Merge everything with sample sizes
tmp2 = left_join(tmp, sz)
# check merge
a1 = sz[grepl("INSEAD", sz$University),]
a2 = dat %>% filter(name=="french")
setdiff(a1$Study, a2$Study)
setdiff(a2$Study, a1$Study)
dat[grepl("Higher Standard", dat$Study),]
dim(tmp2)
apply(tmp2, 2, FUN=function(x) mean(is.na(x)))
dd = tmp2 %>% filter(is.na(Sample.Size)) %>% arrange(University)
sz %>% filter(Study %in% dd$Study)
unique(tmp2$University)
unique(tmp2$Study)


#----Clean up names
out = tmp2 %>% select(experiment=Study, site=name, 
                      n=Sample.Size, d=EffectSize)

#----Compute variances and effect sizes
out$vd = 4/out$n + out$d^2/(2*out$n)
out$g = (1 - 3/(4*(out$n - 2) - 1))*out$d # small sample correction for g
out$vg = (1 - 3/(4*(out$n - 2) - 1))^2*out$vd # small sample correction for g

#---Denote effect size
out$es = "smd"

#---Merge with replication determination from table 1
replicated_df = read.csv('replication_det.csv')
out = left_join(out, replicated_df)

#---Standardize 'replicated' column
out$replicated = as.integer(out$replicated == 'Yes')
head(out)

#---Get original study data
origs = read.csv("./ppir_original.csv")
names(origs)
names(out)
out = rbind(out, origs) %>% arrange(experiment, site)
out$es = 'd'

str(out)

#---Write to file
write.csv(out, "../../ppir.csv", row.names=F)

out %>% group_by(experiment) %>% count()
