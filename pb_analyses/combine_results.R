# set the working directory to src so we can use relative paths
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path), "/results"))

inc = "exclude"

# List of files to load for Q-test of multiple studies
toload = grep("qtest", grep(inc, list.files(), value=T), value=T)



