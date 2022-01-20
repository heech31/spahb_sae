## 1989 Four person median income data analysis using HB spatial model
rm(list=ls())

dataPath <- "./data/"
library(tmvtnorm)
library(boot)
load(paste(dataPath,"dat1989_4.RData",sep="") )

sourcePath <- "./functions/"

# Load packages and source functions
files.sources <-  list.files(path=sourcePath)
extension <- substring(files.sources,nchar(files.sources)-1,nchar(files.sources))
lets.source <- paste(sourcePath, files.sources[extension==".R"], sep="")
mapply(source, lets.source) 


# Scale dollor amount by 1000
dat	<- dat1989_4[,3:6]/1000
head(dat)

# Direct estimates
Y	<- dat[,"y"]
# Sampling errors
sDi	<- dat[,"sDi"]
# Covariates
X	<- as.matrix( cbind(1, dat[,2:3] ) )

# Sampling variances
D <- diag( sDi^2 )
D <- as(D,"sparseMatrix")

# Number of covariates
p <- ncol(X)

# Adjacency matrix
W <- as(W,"sparseMatrix")

rejection.sampling(Y,X,D,W,type="fh",100)
rejection.sampling(Y,X,D,W,type="car",100)
rejection.sampling(Y,X,D,WSAR,type="sar",100)





