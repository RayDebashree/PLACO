#-------- Download to local directory OR directly source PLACO from github
# source("PLACO_v0.1.1.R")
require(devtools)
source_url("https://github.com/RayDebashree/PLACO/blob/master/PLACO_v0.1.1.R?raw=TRUE")
#--------

set.seed(1)
## For an example, let's first simulate a toy set of GWAS summary 
## statistics on 2 uncorrelated traits and 1000 variants
require(MASS)
k <- 2
p <- 1000
Z.matrix <- mvrnorm(n=p, mu=rep(0,k), Sigma=diag(1,k))
P.matrix <- matrix(NA, nrow=p, ncol=k)
for(j in 1:k){
    P.matrix[,j] <- sapply(1:nrow(Z.matrix), 
        function(i) pchisq(Z.matrix[i,j]^2,df=1,ncp=0,lower.tail=F))
}	
colnames(Z.matrix) <- paste("Z",1:k,sep="")
colnames(P.matrix) <- paste("P",1:k,sep="")

## Steps to implementing PLACO
# Step 1: Obtain the variance parameter estimates (only once)
VarZ <- var.placo(Z.matrix, P.matrix, p.threshold=1e-4)
# Step 2: Apply test of pleiotropy for each variant
out <- sapply(1:p, function(i) placo(Z=Z.matrix[i,], VarZ=VarZ))
# Check the output for say variant 100
dim(out)
out[,100]$T.placo
out[,100]$p.placo

## If the traits are dependent or correlated, we suggest
## decorrelating the Z-scores (only once), then apply Steps 1 and 2
## on the decorrelated Z-scores

# Step 0a: Obtain the correlation matrix of Z-scores
R <- cor.pearson(Z.matrix, P.matrix, p.threshold=1e-4)
# Step 0b: Decorrelate the matrix of Z-scores
	# function for raising matrix to any power
	"%^%" <- function(x, pow)
		with(eigen(x), vectors %*% (values^pow * t(vectors)))
Z.matrix.decor <- Z.matrix %*% (R %^% (-0.5))
colnames(Z.matrix.decor) <- paste("Z",1:k,sep="")
