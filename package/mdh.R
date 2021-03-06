mdh_constvar = function(k, alpha=.05, power=0.8, fixed=TRUE, h0replication=TRUE, lambda0=0, tau0=0, step=0.001, maxratio=4){
  
  if(h0replication){
    if(fixed){
      lambda = seq(lambda0/(k-1), max(maxratio, maxratio*lambda0/(k-1)), by=step)*(k-1)
      calpha = qchisq(1-alpha, df=k-1, ncp=lambda0)
      p0 = 1 - pchisq(calpha, df=k-1, ncp=lambda)
      return(min(lambda[which(abs(p0 - power) == min(abs(p0 - power)))])/(k-1))
    } else {
      tau = seq(tau0, maxratio*tau0, by=step)
      
    }
  } else {
    if(lambda0==0){
      print('Error: must specify lambda0 > 0.')
      return(NULL)
    } else {
      lambda = seq(0, lambda0/(k-1), by=step)*(k-1)
      calpha = qchisq(alpha, df=k-1, ncp=lambda0)
      p0 = pchisq(calpha, df=k-1, ncp=lambda)
      return(max(lambda[which(abs(p0 - power) == min(abs(p0 - power)))]/(k-1)))
    }
  }
}


mdh_constvar2 = function(k, alpha=.05, power=0.8, fixed=TRUE, h0replication=TRUE, lambda0=0, tau0=0, epsilon=0.01, maxits=10000000){
  its = 0
  if(h0replication){
    if(fixed){
      rr = 0.5
      lambdai = rr*(k-1)
      calpha = qchisq(1-alpha, df=k-1, ncp=lambda0)
      p0 = 1 - pchisq(calpha, df=k-1, ncp=lambdai)
      while(abs(p0 - power) > epsilon & its < maxits){
        # if(p0 > power){ rr = rr * power/p0 } else { rr = 2*rr }
        rr = power/p0
        lambdai = rr*(k-1)
        calpha = qchisq(1-alpha, df=k-1, ncp=lambda0)
        p0 = 1 - pchisq(calpha, df=k-1, ncp=lambdai)
        its = its + 1
      }
      return(rr)
      
    } else {
      
      return(NULL)
      
    }
  } else {
    
    if(lambda0==0){
      
      print('Error: must specify lambda0 > 0.')
      return(NULL)
      
    } else {
      
      rr = 0.5
      lambdai = 
      calpha = qchisq(alpha, df=k-1, ncp=lambda0)
      p0 = pchisq(calpha, df=k-1, ncp=lambda)
      return(max(lambda[which(abs(p0 - power) == min(abs(p0 - power)))]/(k-1)))
    }
  }
}

# start = Sys.time()
# mdh_constvar(36, maxratio = 100)
# Sys.time() - start
# 
# start = Sys.time()
# mdh_constvar2(36)
# Sys.time() - start


# mdh_constvar(37)
# mdh_constvar(37, lambda0=36/4)
# mdh_constvar(37, lambda0=36/3)
# mdh_constvar(37, lambda0=2/3*36)
# 
# mdh_constvar(37, lambda0=(k-1)/4, BOR=F)
# mdh_constvar(37, lambda0=(k-1)/2, BOR=F)
# mdh_constvar(37, lambda0=(k-1)*2/3, BOR=F)
# mdh_constvar(37, lambda0=(k-1), BOR=F)
