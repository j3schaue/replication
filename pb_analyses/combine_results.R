# set the working directory to src so we can use relative paths
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path), "/results"))

inc = "exclude"

# List of files to load for Q-test of multiple studies
toload = c(grep("qtest", grep(inc, list.files(), value=T), value=T),
           "qtest_fixed_rpe_include.csv", "qtest_fixed_rpp_include.csv")
dfs = lapply(toload, read.csv)

res = do.call(rbind, dfs)

res %>% group_by(paper) %>% 
  summarize(total = n(), nonrep0 = sum(p0 < 0.05), nonrep25 = sum(p25 < 0.05), 
            nonreplication = sum(replicated == 0), 
            mdh0 = mean(mdh0), mdh25 = mean(mdh25))
