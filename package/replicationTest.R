###--------------------------------------------###
###--------------------------------------------###
### TESTING REPLICATION
###
### This code provides some useful functions for
### testing whether results from a series of 
### experiments confirm or disconfirm replication.
###--------------------------------------------###
###--------------------------------------------###

library(metafor)

###--------------------------------------------###
### Calculating the test statistic
###
### The Q-statistic is used for most of the tests.
### We assume that findings are normally 
### distributed around an unknown, study-specific
### mean theta_i, with known sampling variance v_i.
###--------------------------------------------###

weightedMean <- function(data=NULL, t=NULL, v=NULL, random=F){
  # Takes: data (data.frame) with two columns 't', and 'v'
  #         OR t (vector of effect sizes) and v (vector of variances)
  #
  # Returns: precision-weighted mean of effect sizes (numeric)
  
  if(!is.null(data)){ # if the effect sizes and variances are passed as a data.frame, set them as vectors.
    t <- as.numeric(data$t)
    v <- as.numeric(data$v)
  }
  
  t_bar_dot <- sum(t/v)/sum(1/v)
  return(t_bar_dot)
}

qStat <- function(data=NULL, t=NULL, v=NULL){
  # Takes: data (data.frame) with two columns 't', and 'v'
  #         OR t (vector of effect sizes) and v (vector of variances)
  #
  # Returns: Q statistic (numeric)
  
  t_bar_dot <- weightedMean(data, t, v)
  
  if(!is.null(data)){
    t <- data$t
    v <- data$v
  }
  
  Q <- sum((t - t_bar_dot)^2 / v)
  return(Q)
}


###--------------------------------------------###
### Hypothesis test 
###--------------------------------------------###

replicationTest <- function(data=NULL, t=NULL, v=NULL, burden_on_replication=TRUE, fixed=TRUE, exact=TRUE, alpha=.05, lambda0=NULL, tau0=NULL){
  
  Q <- qStat(data, t, v) # Q statistic
  k <- max(nrow(data), length(t)) # number of studies
  g <- NULL
  h <- NULL
  
  if(!is.null(data)){
    t <- data$t
    v <- data$v
  }
  
  if(exact){ # Exact replications
            ## 2 Cases: approximate, fixed, nonreplication burden; approximate, randome, nonreplication burden
    if(burden_on_replication){ # No test if the burden is on replication (i.e., null is that studies do not replicate)
      print("There is no test for exact repliciation under this null hypothesis.")
      return()
    } else { # If null is that studies replicate, then  Q has the same null distribution whether or not we consider stuides as fixed.
      calpha <- qchisq(1 - alpha, df=k-1, ncp=0)
      reject <- Q >= calpha
      p <- 1 - pchisq(Q, df=k-1, ncp=0)
    }
  } else if(fixed & !exact){ ### CASES: Studies are fixed, replication is approximate
      if(is.null(lambda0)){ # If there is no pre-specied heterogeneity parameter, use (k-1)/3 (Hunter & Schmidt)
        lamdbda0 <- (k - 1)/3
        print("No parameter relating to feasible heterogeneity specified. Using lambda=(k-1)/3.")
      }
      if(burden_on_replication){
        calpha <- qchisq(alpha, df=k-1, ncp=lambda0)
        reject <- Q <= calpha
        p <- pchisq(Q, df=k-1, ncp=lambda0)
      } else {
        calpha <- qchisq(1 - alpha, df=k-1, ncp=lambda0)
        reject <- Q >= calpha
        p <- 1 - pchisq(Q, df=k-1, ncp=lambda0)
      }
  } else if (!fixed & !exact) {### CASES: Studies are random, replication is approximate
      if(is.null(tau0)){ # If there is no pre-specified hetereogeneity, use v/2
        tau0 <- 0.5 * mean(v)
      }
      
      # Parameters for the Gamma approximation (Hedges & Pigott; Satterthwaite)    
      const <- sum(1/v) - sum(1/v^2)/sum(1/v)
      muq <- k - 1 + const * tau0
      varq <- 2*(k - 1) + 4 * const * tau0 + 2*(sum(1/v^2) - 2*sum(1/v^3)/sum(1/v) + sum(1/v^2)^2/sum(1/v)^2)*tau0^2
      g <- varq/(2*muq)
      h <- 2*muq^2/varq
      
      if(burden_on_replication){
        calpha <- varq/(2*muq) * qchisq(alpha, df=2*muq^2/varq, ncp=0)
        reject <- Q <= calpha
        p <- pchisq(Q/g, df=h, ncp=0)
      } else {
        calpha <- varq/(2*muq) * qchisq(1 - alpha, df=2*muq^2/varq, ncp=0)
        reject <- Q >= calpha
        p <- 1 - pchisq(Q/g, df=h, ncp=0)
      }
  }
  
#   print(cat(
#         paste0("The test statistic is Q=", round(Q, digits=3)),
#         paste0("The critical value is ", round(calpha, digits=3)),
#         paste0("The test ", ifelse(reject, "rejects", "fails to reject"), " the null hypothesis (p=", round(p, digits=3), ")."),
#         paste0("-----------------------------"),
#         paste0("Data: ", k, " studies"),
#         paste0("Null hypothesis: ", ifelse(burden_on_replication, "Studies do not replicate", "Studies replicate")),
#         paste0("Studies treated as ", ifelse(fixed, "fixed.", "random.")),
#         paste0("Testing for ", ifelse(exact, "exact ", "approximate "), "replication."), 
#         sep="\n"
#         ))
  
  return(list(stat=Q, p=p, reject=reject, crit=calpha, h=h, g=g))
}

# k <- 30
# v <- runif(k, .1, .2)
# mn <- rnorm(k, .2, .351)
# data <- data.frame(t=rnorm(k, mn, v), v=v)
# replicationTest(data, burden_on_replication=T, fixed=F, exact=F, tau0=mean(v)/2)


###--------------------------------------------###
###--------------------------------------------###
### Power Functions
###--------------------------------------------###
###--------------------------------------------###
# Note that lambda0 for approximate replications is expressed 
# in terms of the ratio of between to within variance, (k-1) tau/v
# For tau0 for approximate replication, we express it as a 
# proportion of v.

###-----Burden of proof is on nonreplication

pwrNRFE <- function(k, lambda, alpha=.05){
  # burden on non-replication, fixed studies, exact replication
  
  c_alpha <- qchisq(1 - alpha, df=k-1, ncp=0)
  return(1 - pchisq(c_alpha, df=k-1, ncp=lambda))
}

pwrNRFA <- function(k, lambda, lambda0=NULL, alpha=.05){
  # burden on non-replication, fixed studies, approximate replication
  
  if(is.null(lambda0)){# if no null value provided, choose it to be half of the true value
    lambda0 <- lambda/2
  }
  
  c_alpha <- qchisq(1 - alpha, df=k-1, ncp=lambda0)
  return(1 - pchisq(c_alpha, df=k-1, ncp=lambda))
}

pwrNRRE <- function(k, tau, v=1, alpha=.05){
  # burden on non-replication, random studies, exact replication
  
  c = ifelse(length(v)==1, (k-1)/v, sum(1/v) - sum(1/v^2)/sum(1/v))
  d = ifelse(length(v==1), (k-1)/v^2, sum(1/v^2) - 2*sum(1/v^3)/sum(1/v) + sum(1/v^2)^2/sum(1/v)^2)
  muq <- (k - 1)  + c*tau
  varq <- 2*(k - 1)  + 4*c*tau + 2*d*tau^2
  
  g <- varq/(2 * muq)
  h <- 2 * muq^2 / varq

  c_alpha <- qchisq(1 - alpha, df=k-1, ncp=0)
  return(1 - pchisq(c_alpha/g, df=h, ncp=0))
}

pwrNRRA <- function(k, tau, tau0=NULL, v=1, alpha=.05){
  # burden on non-replication, random studies, approximate replication
  
  if(is.null(tau0)){ # if no null value tau0 is entered, assume it is tau/2
    tau0 <- tau/2
  }
  
  c = ifelse(length(v)==1, (k-1)/v, sum(1/v) - sum(1/v^2)/sum(1/v))
  d = ifelse(length(v==1), (k-1)/v^2, sum(1/v^2) - 2*sum(1/v^3)/sum(1/v) + sum(1/v^2)^2/sum(1/v)^2)
  
  muq0 <- (k - 1)  + c*tau0
  varq0 <- 2 * (k - 1)  + 4*c*tau0 + 2*d*tau0^2
  g0 <- varq0/(2 * muq0)
  h0 <- 2 * muq0^2 / varq0
  
  muq <- (k - 1)  + c*tau
  varq <- 2 * (k - 1)  + 4*c*tau + 2*d*tau^2
  g <- varq/(2 * muq)
  h <- 2 * muq^2 / varq
  
  c_alpha <- g0*qchisq(1 - alpha, df=h0, ncp=0)
  return(1 - pchisq(c_alpha/g, df=h, ncp=0))
}

pwrRFA <- function(k, lambda, lambda0, v=1, alpha=.05){
  # burden on replication, fixed studies, approximate replication
  
  c_alpha <- qchisq(alpha, df=k-1, ncp=lambda0)
  return(pchisq(c_alpha, df=k-1, ncp=lambda))
}

pwrRRA <- function(k, tau, tau0, v=1, alpha=.05){
  # burden on replication, random studies, approximate replication
  
  c = ifelse(length(v)==1, (k-1)/v, sum(1/v) - sum(1/v^2)/sum(1/v))
  d = ifelse(length(v==1), (k-1)/v^2, sum(1/v^2) - 2*sum(1/v^3)/sum(1/v) + sum(1/v^2)^2/sum(1/v)^2)
  
  muq0 <- (k - 1)  + c*tau0
  varq0 <- 2 * (k - 1)  + 4*c*tau0 + 2*d*tau0^2
  g0 <- varq0/(2 * muq0)
  h0 <- 2 * muq0^2 / varq0
  
  muq <- (k - 1)  + c*tau
  varq <- 2 * (k - 1)  + 4*c*tau + 2*d*tau^2
  g <- varq/(2 * muq)
  h <- 2 * muq^2 / varq
  
  calpha <- g0*qchisq(alpha, df=h0, ncp=0)
  return(pchisq(calpha/g, df=h, ncp=0))
}

