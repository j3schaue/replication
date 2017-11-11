# Function for converting a column of Qualtrics output to useable form #
convert <- function(x) as.numeric(as.character(x)) 

# Function for Data Summary with NA's 
my_summary <- function(v){
  if(!any(is.na(v))){
    res <- c(summary(v),"NA's"=0)
  } else res <- summary(v)
  return(res)
}

# Raw Data Summary for Each Lab #
csv <- list.files(pattern="*.csv") # names of all csv files in folder
authors <- sapply(strsplit(csv,split='_',fixed=TRUE),function(x){x[1]}) # authors' names
k <- length(authors) # number of labs
data.all <- list()
vars <- c(18,22,24,28,31,32,35:78,81:86,88)
for(i in 1:k){
  data <- read.csv(csv[i],header=FALSE) # read csv file
  if(as.character(data[1,1])=='startDate'){ # if data in new format
    data <- data[,c(1:10,14:22,11,23:27,12,28:30,13,31:33,13,34:87)][,vars] # make compatible
  } else data <- data[,vars] # otherwise, no change
  colnames(data) <- sapply(data[1,],as.character) # label columns
  data <- suppressWarnings(sapply(data[-(1:2),],convert))
  data[,c(2,4)][data[,c(2,4)]==0] <- NA
  sum.dat <- t(apply(data,2,my_summary))
  write.csv(sum.dat,paste(authors[i],'_Summary.csv',sep=''))
}
