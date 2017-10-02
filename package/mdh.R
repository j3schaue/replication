mdh_constvar = function(k, alpha=.05, power=0.8, BON=TRUE, lambda0=0, step=0.001, maxratio=4){
  
  if(BON){
    lambda = seq(lambda0/(k-1), max(maxratio, maxratio*lambda0/(k-1)), by=step)*(k-1)
    calpha = qchisq(1-alpha, df=k-1, ncp=lambda0)
    p0 = 1 - pchisq(calpha, df=k-1, ncp=lambda)
    # print(calpha); print(k); print(lambda0)
    return(min(lambda[which(abs(p0 - power) == min(abs(p0 - power)))])/(k-1))
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

# mdh_constvar(37)
# mdh_constvar(37, lambda0=36/4)
# mdh_constvar(37, lambda0=36/3)
# mdh_constvar(37, lambda0=2/3*36)
# 
# mdh_constvar(37, lambda0=(k-1)/4, BOR=F)
# mdh_constvar(37, lambda0=(k-1)/2, BOR=F)
# mdh_constvar(37, lambda0=(k-1)*2/3, BOR=F)
# mdh_constvar(37, lambda0=(k-1), BOR=F)
