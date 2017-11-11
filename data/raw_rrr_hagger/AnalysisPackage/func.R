# Setup a useful function for calculating response when there are two ePrime objects where response may occur
msit.RespCalc=function(stimRESP,jitRESP,srcind){
  resp=NA
  src=''
  flag=0
  stimRESP=as.numeric(stimRESP)
  jitRESP=as.numeric(jitRESP)
  if (!is.na(stimRESP)){resp=stimRESP;src='stim'} else if(!is.na(jitRESP)) {resp=jitRESP;src='jit'}
  if(srcind=='resp'){return(resp)}
  if(srcind=='src'){return(src)}
}

LetE.RespCalc=function(stimRESP,jitRESP,srcind){
  resp=0
  src=''
  stimRESP=as.numeric(stimRESP)
  jitRESP=as.numeric(jitRESP)
  if (stimRESP==1){resp=stimRESP;src='stim' # if positSpectrum and spectral density estimation by the Discrete Fourier transform (DFT), including a comprehensive list of window functions and some new at-top windowsive response in stim, use it. if not, go on
  } else if(jitRESP==1) # is positive response in jitterfix, use it. if not, stop, and use default 0 value
  {resp=jitRESP;src='jit'}
  if(srcind=='resp'){return(resp)}
  if(srcind=='src'){return(src)}
}

int.trapz <- function(y,x=(1:length(y))){
  idx <- 2:length(x)
  val <-  ( (x[idx] - x[idx-1]) %*% (y[idx] + y[idx-1]) ) / 2
  return(val)
}

win.daniell <- function(N){
  win=c()
  for (i in 1:N){
    if(i==1 | i==N){
      win[i]=1/(2*(N-1))
    } else {
      win[i]=1/(N-1)
    }
  }
  return(win)
}

win.hanning <- function(N){
  win <- (.5 * (1-cos(2*pi*(0:(N-1))/(N-1))))
  return(win/sum(win))
}


         
WindowApply <- function(data,WinSize,func,...){
  # input arguments
  ## data - dataframe
  ## WinSize - How big the window should be
  ## func - function handle
  ## ... - add'l arguments to func
  # Operation
  ## Operate row-wise on dataframe
  ## Loop rowwise, split row vector into moving windows, apply func

nWin <- ncol(data) - WinSize + 1
out = data.frame(matrix(nrow=nrow(data),ncol=ncol(data)))


  for (iR in 1:nrow(data)){
    curData = as.numeric(data[iR,])
    for (iW in 1:nWin){
      out[iR,iW] = func(curData[iW:(iW+WinSize-1)],...)
    }
  }
return(out)
}
    
mexgauss2 <- function(x,arg){
  x = x[!is.na(x)]
  if ( length(x) > 1){
    ExG <- mexgauss(x)
    return(ExG[arg])
  } else return(NA)
}
  
