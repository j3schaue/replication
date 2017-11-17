
# set the working directory to src so we can use relative paths
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

dfs = list(
  rpp = read.csv("./results/qtest_fixed_rpp.csv"),
  rpe = read.csv("./results/qtest_fixed_rpe.csv"),
  manlabsf = read.csv("./results/comparison_manylabs_FE.csv"),
  manylabsr = read.csv("./results/comparison_manylabs_DL.csv"),
  alognaf = read.csv("./results/comparison_alogna_FE.csv"),
  alognar = read.csv("./results/comparison_alogna_DL.csv")
)

lambda0s = c(0, 1/4, 1/3, 2/3)
strlambdas = round(100*lambda0s, 0)

prop_replicated = setNames(data.frame(matrix(unlist(lapply(strlambdas, FUN=function(pp)
  sapply(dfs, FUN=function(df) mean(df[[paste0("p", pp)]] > 0.05, na.rm=T))
)), ncol = length(dfs), byrow=T)), names(dfs))

prop_replicated$lambda0 = lambda0s

agree = setNames(
  data.frame(
    matrix(unlist(
          lapply(strlambdas, FUN=function(pp)
            sapply(dfs, FUN=function(df){
              mean(df$replicated == as.integer(df[[paste0('p', pp)]] > .05), na.rm=T)
          })
)), ncol=length(dfs), byrow=T)), names(dfs))

agree$lambda0 = lambda0s
agree

# original studies say that findings don't replicate, 
# but the Q test is inconclusive
falsepos = setNames(
  data.frame(
    matrix(unlist(
      lapply(strlambdas, FUN=function(pp)
        sapply(dfs, FUN=function(df){
          mean((df$replicated == 0 & df[[paste0('p', pp)]] > .05), na.rm=T)
        })
      )), ncol=length(dfs), byrow=T)), names(dfs))

falsepos$lambda0 = lambda0s
falsepos

# Original studies say that studies replicate, but
# the Q test disagrees
falseneg = setNames(
  data.frame(
    matrix(unlist(
      lapply(strlambdas, FUN=function(pp)
        sapply(dfs, FUN=function(df){
          mean((df$replicated == 1 & df[[paste0('p', pp)]] < .05), na.rm=T)
        })
      )), ncol=length(dfs), byrow=T)), names(dfs))

falseneg$lambda0 = lambda0s
falseneg

for(pp in strlambdas){
  assign(paste0("disagrees", pp),
         lapply(dfs, FUN=function(df){
            df[which(df$replicated == as.integer(df[[paste0("p", pp)]] < 0.05)), ]
         })
  )
}
  
for(i in 1:length(disagrees0)){
  print(ncol(disagrees0[[i]]))
  if(nrow(disagrees0[[i]]) > 0){
    disagrees0[[i]]$paper = names(disagrees0)[i] 
  }
}
disagrees0$rpe[, c("p0", "replicated")]
