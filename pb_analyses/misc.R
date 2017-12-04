
combineResults = function(t=NULL, v=NULL, h0replication=TRUE, fixed=TRUE, alpha=.05, lambda0=0, tau0=0, power=0.8, step=.001, maxratio=100){
  qtest = replicationTest(t=t, v=v, h0replication=h0replication, fixed=fixed, alpha=alpha, lambda0=lambda0, tau0=tau0)
  qtest[["mdh"]] = mdh_constvar(k=length(t), alpha=alpha, power=power, h0replication=h0replication, lambda0=lambda0, step=step, maxratio=maxratio)
  return(qtest[c("k", "Q", "calpha", "p", "mdh")])
}

runComparisonAnalyses = function(data, t, v, ratios, paper, methods){
  require(metafor)
  
  # List of experiments
  experiments = unique(data$experiment)
  
  for(mm in methods){# loop through methods
    
    # comp is a table for each experiment and lambda0 for a given synthetic replicate.
    comp = lapply(experiments, FUN=function(ee){
      print(ee)
      
      # separate out the replicates and original study
      replicates = data %>% filter(experiment==ee & site!='original')
      orig = data %>% filter(experiment==ee & site=='original')
      
      # combine the replicates
      tmp = rma.uni(yi=replicates[[t]], vi=replicates[[v]], method=mm)
      
      foo = lapply(ratios, FUN=function(rr){ # loop through the lambda0s
        
        lambda0 = rr # set the null hypothesis
        
        if(nrow(orig) == 1){
          
          # run a Q-test and get the MDH
          ll = data.frame(combineResults(t=c(tmp$beta[,1], orig[[t]]), 
                                         v=c(tmp$se^2, orig[[v]]),
                                         lambda0=lambda0, maxratio=20))
          
        } else {
          ll = data.frame(k = 2, Q = NA, calpha = NA, p = NA, mdh = NA)
        }
        
        for(i in 3:5){
          names(ll)[i] = paste0(names(ll)[i], round(rr*100, 0))
        }
        
        return(ll)
      })
      
      bar = Reduce(left_join, foo)
      if(nrow(orig) == 0) {
        bar$t1 = NA; bar$v1 = NA
      } else {
        bar$t1 = orig[[t]]; bar$v1 = orig[[v]]
      }
      bar$t2 = tmp$beta[1,1]; bar$v2 = tmp$se^2
      
      return(bar)
    })
    
    # aggregate results
    comptab = Reduce(rbind, comp)
    comptab$experiment = experiments
    comptab$paper = paper
    
    print(names(comptab))
    
    # clean up order of columns and write to file
    comptab = comptab[c("paper", "experiment", "k", "Q", 
                        "calpha0",  "p0", "mdh0", 
                        "calpha25", "p25",  "mdh25", 
                        "calpha33", "p33",  "mdh33", 
                        "calpha67", "p67", "mdh67", 
                        "t1", "t2", "v1", "v2")] %>%
      left_join(., select(data, experiment, replicated)) %>% distinct()
    
    write.csv(comptab, paste0("./results/comparison_", paper, "_", mm, ".csv"), row.names=F)
  }
}


runQtests = function(data, ratios, t, v, paper, include=FALSE){
  
  ## Set parameters for analysis
  experiments = unique(data$experiment) # unique experiment names
  if('n' %in% names(data)){
    ks = sapply(experiments, # no. of trials per experiment
                FUN=function(ee) count(filter(data, experiment==ee))$n)
  } else {
    ks = sapply(experiments, # no. of trials per experiment
                FUN=function(ee) count(select(filter(data, experiment==ee), -n))$n)
  }
  vbars = sapply(experiments, FUN=function(ee) mean(dplyr::filter(data, experiment==ee)$vd)) # avg sampling variances
  
  fe = lapply(ratios, FUN=function(rr) # loop through null hypotheses 
    setNames(data.frame( # store results as a data frame
      matrix(unlist(lapply(seq_along(experiments), FUN=function(i)
        
        # get the results of the Q-test using all of the studies (rather than aggregating replicates)
        combineResults(t=filter(data, experiment==experiments[i])[[t]],
                       v=filter(data, experiment==experiments[i])[[v]],
                       lambda0=(ks[i]-1)*rr, 
                       maxratio=20)
      )
      ), ncol=5, byrow = T)
    ), c("k", "Q", paste0("calpha", round(rr*100, 0)), # name the columns
         paste0("p", round(rr*100, 0)), paste0("mdh", round(rr*100, 0)))
    )
  )
  
  # Join for all null hypotheses
  fetab = Reduce(left_join, fe)
  fetab$experiment = experiments
  fetab$vbar = vbars
  fetab$paper = paper
  
  # reorder data and write to file
  fetab = fetab[c("paper", "experiment", "k", "Q", "calpha0",  "p0", "mdh0", "calpha25", "p25",  "mdh25", "calpha33", "p33",  "mdh33", "calpha67", "p67", 
                  "mdh67", "vbar")] %>%
    left_join(., distinct(dplyr::select(data, experiment, replicated)))
  
  write.csv(fetab, "./results/qtest_fixed_", paper, "_include.csv", row.names=F)
}



qCI = function(Q, k, alpha=.05){
  lb = 0
  pct = 1 - pchisq(Q, k-1)
  if(pct <  alpha){
    lb = optim(par=0, fn = function(x) abs(pchisq(Q, k-1, x) - 1 + alpha))$par
  }
  
  ub = Q
  pct = pchisq(Q, k-1, ub)
  if(pct < alpha){
    ub = optim(par=0, fn = function(x) abs(pchisq(Q, k-1, x) - alpha))$par
  }
  
  return(list(lblambda=lb, ublambda=ub))
}
