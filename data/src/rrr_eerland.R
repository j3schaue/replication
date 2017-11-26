library(dplyr)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Loading raw data
im <- read.csv("../raw/rrr_eerland/Forest plots R/imagery.csv")
im["experiment"] <- "Imagery"
im["es"] <- "d"
im["replicated"] <- "0"
im <- rename(im, site = Author)

ia <- read.csv("../raw/rrr_eerland/Forest plots R/Intention_attribution.csv")
ia["experiment"] <- "Intention_attribution"
ia["es"] <- "d"
ia["replicated"] <- "0"
ia <- rename(ia, site = Author)

intent <- read.csv("../raw/rrr_eerland/Forest plots R/Intentionality.csv")
intent["experiment"] <- "Intentionality"
intent["es"] <- "d"
intent["replicated"] <- "0"
intent <- rename(intent, site = Author)

#Generating Meta-Analysis Effect Sizes
im <- mutate(im,
             d = (mimp-mperf)/sqrt(((nimp-1)*(sdimp)^2+(nperf-1)*(sdperf)^2)/(nimp+nperf-2)),
             vd = (nimp+nperf)/(nimp*nperf)+(d^2)/(2*(nimp+nperf)),
             df = nimp+nperf-2)

ia <- mutate(ia,
             d = (mimp-mperf)/sqrt(((nimp-1)*(sdimp)^2+(nperf-1)*(sdperf)^2)/(nimp+nperf-2)),
             vd = (nimp+nperf)/(nimp*nperf)+(d^2)/(2*(nimp+nperf)),
             df = nimp+nperf-2)
intent <- mutate(intent,
             d = (mimp-mperf)/sqrt(((nimp-1)*(sdimp)^2+(nperf-1)*(sdperf)^2)/(nimp+nperf-2)),
             vd = (nimp+nperf)/(nimp*nperf)+(d^2)/(2*(nimp+nperf)),
             df = nimp+nperf-2)

#Combining experiments
df <- bind_rows(im, ia, intent)
df$site <- gsub('H & A EXPERIMENT 3', 'original', df$site)
df$site[grepl('ONLINE', df$site)] = 'online'

#Writing CSV file
write.csv(df, "../rrr_eerland.csv", row.names=F)
