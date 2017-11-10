###-----------------------------------------------------------###
###-----------------------------------------------------------###
### PPIR files
###-----------------------------------------------------------###
###-----------------------------------------------------------###

library(dplyr); library(tidyr)

setwd("~/Documents/replication/replication/data/raw_ppir/") #path on Jake's machine
# setwd("./data/raw_ppir") # relative path


sz = read.csv("schweinsberg_samplesizes_new.csv")
names(sz)[1] = "Study"
sz$Study = as.character(sz$Study)
sz$University = as.character(sz$University)


sz$University = gsub(" in Washington DC", "",
  gsub(", (Canada|Germany|the Netherlands|China)", "", sz$University))

nm = read.csv("name-uni.csv")
nm$name = as.character(nm$name)
nm$uni = as.character(nm$uni)
nm$name = gsub("jenniferjordan", "jenjordan",
            gsub("victoriabrescoll", "toribrescoll",
              gsub("xiaominsun", "sunnysun",
                gsub("danmolden[[:alnum:]]+", "danmolden",
                  tolower(gsub("[[:blank:]]", "", nm$name))))))
nmsun = sort(unique(nm$name))
nm[nm$name=="danmolden",]$uni = "Northwestern University"
nm = unique(nm)
nm[nm$name=="inseadsorbonnelab",] = c("french", "INSEAD, France")
nm = rbind(nm, data.frame(name=c("mturk", "ucibusiness", "ucipsychstudents"),
           uni=c("Mechanical Turk sample", "University of California Irvine", "University of California Irvine"))
)

dat = read.csv("Shweinsberg_fig1.csv")
dat$Study = as.character(dat$Study)
dat$Study[dat$Study=="Intuitive economics"] = "Intuitive Economics"
dat$Study[dat$Study=="Cold Hearted Prosociality"] = "Cold-Hearted Prosociality"
dat$name = 
  gsub("annlaure", "annelauresellier",
  gsub("toribrescolli", "toribrescoll", 
   gsub("additionalfrench", "french",
     gsub("data|translation", "",
       tolower(gsub("[[:blank:]]", "", dat$Sample))))))

dat %>%
  filter(name %in% c("french", "frenchtranslation", "additionalfrenchdata"))

sz[grepl("INSEAD", sz$University), ]

setdiff(unique(nm$name), unique(dat$name))
setdiff(unique(dat$name), unique(nm$name))


uni1 = sort(unique(sz$University))
uni2 = sort(unique(nm$uni))
setdiff(uni1, uni2)
setdiff(uni2, uni1)

filter(dat, name=="felixcheung")
filter(sz, University=="University of Hong Kong")

filter(sz, University=="University of Michigan")
filter(dat, name=="tatianasokolova")
filter(sz, University=="HEC Paris")
filter(dat, name=="annelauresellier")

# Merge dat and sz, but we need to connect the dat to universities 
names(nm)[2] = "University"
tmp = left_join(dat, nm) #%>% left_join(., dat)
tmp %>% group_by(Study, University) %>% tally() %>% filter(n>1)
tmp %>% filter(University == "University of Washington (Foster)")
apply(tmp, 2, FUN=function(x) mean(is.na(x)))

dim(tmp)
head(tmp)

tmp2 = left_join(tmp, sz)
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
# tmp %>% filter(University != "University of California Irvine") %>%
#   group_by(Study, name, University) %>% 
#   tally() %>% nrow()
#   arrange(Study, University)
# tmp %>% filter(University != "University of California Irvine") %>%
#   group_by(Study, University) %>% tally() %>% nrow()
#   arrange(Study, University)
# tmp %>% group_by(Study, name) %>% tally() %>% nrow()

head(tmp2)
write.csv(tmp2, "../ppir_full.csv", row.names=F)


out = tmp2 %>% select(experiment=Study, lab=name, 
                      d=EffectSize, n=Sample.Size)
out$vd = 4/out$n + out$d^2/(2*out$n)

write.csv(out, "../ppir.csv", row.names=F)



pp1 = read.csv("PIPR 1.csv")
