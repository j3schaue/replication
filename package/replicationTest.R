###--------------------------------------------------------###
###--------------------------------------------------------###
### TESTING REPLICATION
###
### This code provides some useful functions for
### testing whether results from a series of 
### experiments confirm or disconfirm replication.
###--------------------------------------------------------###
###--------------------------------------------------------###

###--------------------------------------------###
### Calculate the Q-statistic 
###--------------------------------------------###
qStat = function(t=NULL, v=NULL, data=NULL){
  ## Takes: data (data.frame) with two columns 't', and 'v'
  ##         OR t (vector of effect sizes) and v (vector of variances)
  ##
  ## Returns: Q statistic (numeric)
  
  #---sanitize inputs
  if(!is.null(data)){
    t <- data$t
    v <- data$v
  }
  
  # weighted mean
  t_bar_dot = weighted.mean(x=t, w=1/v/sum(1/v))
  
  # return Q
  return(sum((t - t_bar_dot)^2 / v))
}


###--------------------------------------------###
### Hypothesis test 
###--------------------------------------------###
replicationTest = function(t=NULL, v=NULL, data=NULL, h0replication=TRUE, fixed=TRUE, alpha=.05, lambda0=0, tau0=0){
  ## Takes: t, effect sizes, vector
  ##        v, variances, vector
  ##        data, data frame containing columns t and v
  ##        h0replication, the null hypothesis, boolean
  ##           TRUE = h0 is replication, FALSE = h0 is nonreplication 
  ##        fixed, fixed vs random effects, boolean
  ##        alpha, level of the test, numeric between 0 and 1
  ##        lambda0, null hypothesis value for fixed effects, 
  ##            default is to 0 for exact replication
  ##        tau0, null hypothesis value for random effects
  ## Returns:

  #---sanitize input
  if(!is.null(data)){
    t <- data$t
    v <- data$v
  }
  
  Q <- qStat(t, v) # Q statistic
  k <- length(t) # number of studies
  
  if(fixed){ # Fixed effects test
    
    if(h0replication){ # H0: lambda <= lambda0
      calpha = qchisq(1-alpha, df=k-1, ncp=lambda0) # critical value
      reject = Q >= calpha # do we reject H0?
      pval = 1 - pchisq(Q, df=k-1, ncp=lambda0) # p-value
    } else { # H0: lambda >= lambda0
      calpha = qchisq(alpha, df=k-1, ncp=lambda0) # critical value
      reject = Q <= calpha # do we reject H0?
      pval = pchisq(Q, df=k-1, ncp=lambda0) # p-value
    }
    
    out = list(k=k, Q=Q, calpha=calpha, p=pval, alpha=alpha, reject=reject)
  
  } else { # Random effect test
    
    # Parameters for the Gamma approximation (Hedges & Pigott; Satterthwaite)    
    const <- sum(1/v) - sum(1/v^2)/sum(1/v)
    muq <- k - 1 + const * tau0
    varq <- 2*(k - 1) + 4 * const * tau0 + 2*(sum(1/v^2) - 2*sum(1/v^3)/sum(1/v) + sum(1/v^2)^2/sum(1/v)^2)*tau0^2
    g <- varq/(2*muq)
    h <- 2*muq^2/varq
    
    if(h0replication){ # H0: tau <= tau0
      calpha <- varq/(2*muq) * qchisq(1 - alpha, df=2*muq^2/varq, ncp=0) # critical value
      reject <- Q >= calpha # do we reject H0?
      pval <- 1 - pchisq(Q/g, df=h, ncp=0) # p-value
    } else { # H0: tau >= tau0
      calpha <- varq/(2*muq) * qchisq(alpha, df=2*muq^2/varq, ncp=0) # critical value
      reject <- Q >= calpha # do we reject H0?
      pval <- pchisq(Q/g, df=h, ncp=0) # p-value
    }
    
    out = list(k=k, Q=Q, calpha=calpha, p=pval, g=g, h=h, alpha=alpha, reject=reject)
  }
  
  print(cat(
        paste0("The test statistic is Q=", round(Q, digits=3)),
        paste0("The critical value is ", round(calpha, digits=3)),
        paste0("The test ", ifelse(reject, "rejects", "fails to reject"), " the null hypothesis (p=", round(pval, digits=3), ")."),
        paste0("-----------------------------"),
        paste0("Data: ", k, " studies"),
        paste0("Null hypothesis: ", ifelse(h0replication, "Studies replicate", "Studies do not replicate")),
        paste0("Studies treated as ", ifelse(fixed, "fixed.", "random.")),
        paste0("Testing for ", ifelse(tau0==0 & lambda0==0, "exact ", "approximate "), "replication."),
        sep="\n"
        ))
  
  return(out)
}

# k <- 30
# v <- runif(k, .1, .2)
# mn <- rnorm(k, .2, .351)
# data <- data.frame(t=rnorm(k, mn, v), v=v)
# replicationTest(data=data, h0replication=T, fixed=F, tau0=mean(v)/2)

###--------------------------------------------###
### Power
###--------------------------------------------###
powerRepTest = function(k, v, fixed=TRUE, h0replication=TRUE, alpha=.05, lambda=NULL, lambda0=0, tau=NULL, tau0=0){
  
  if(fixed){
    if(h0replication){
      pwr = 1 - pchisq(qchisq(1-alpha, df=k-1, ncp=lambda0), df=k-1, ncp=lambda)
    } else{
      pwr = pchisq(qchisq(alpha, df=k-1, ncp=lambda0), df=k-1, ncp=lambda)
    }
  } else {
    cnst = sum(1/v) - sum(1/v^2)/sum(1/v)
    denom = sum(1/v^2) - 2*sum(1/v^3)/sum(1/v) + sum(1/v^2)^2/sum(1/v)^2
    
    # parameters for null distribution
    muq0 <- (k - 1)  + cnst*tau0
    varq0 <- 2 * (k - 1)  + 4*cnst*tau0 + 2*denom*tau0^2
    g0 <- varq0/(2 * muq0)
    h0 <- 2 * muq0^2 / varq0
    
    # parameters for non-null distribution
    muq <- (k - 1)  + cnst*tau
    varq <- 2 * (k - 1)  + 4*cnst*tau + 2*denom*tau^2
    g <- varq/(2 * muq)
    h <- 2 * muq^2 / varq
    
    if(h0replication){
      pwr = 1 - pchisq(g0*qchisq(1 - alpha, df=h0, ncp=0)/g, df=h)
    } else {
      pwr = pchisq(g0*qchisq(alpha, df=h0, ncp=0)/g, df=h)
    }
  }
  return(pwr)
}
# powerRepTest(k=20, lambda=2/3) # exact, fixed
# powerRepTest(k=20, lambda=2/3, lambda0=1/4) # approximate, fixed
# powerRepTest(k=20, h0replication=FALSE, lambda=0, lambda0=1/4) # approximate, fixed, equivalence


###--------------------------------------------###
### MDH/MAH
###--------------------------------------------###
