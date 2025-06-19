#-------- Download to local directory OR directly source PLACO+ from github
# source("PLACO_v0.2.0.R")
require(devtools)
source_url("https://github.com/RayDebashree/PLACO/blob/master/PLACO_v0.2.0.R?raw=TRUE")
#--------

set.seed(1)
## For an example, let's first simulate a toy set of GWAS summary 
## statistics on 2 correlated traits and 1000 variants

# Load the necessary library
if (!require(MASS)) install.packages("MASS")
library(MASS)

# Parameters for simulation
k <- 2    # Number of traits
p <- 1000 # Number of variants
r <- 0.5  # Desired correlation between traits

# Covariance matrix with specified correlation
Sigma <- matrix(c(1, r, r, 1), ncol = 2)

# Generate Z-scores for two correlated traits across 1000 variants
Z.matrix <- mvrnorm(n=p, mu=rep(0,k), Sigma=Sigma)

# Calculate p-values from Z-scores
P.matrix <- matrix(NA, nrow=p, ncol=k)
for(j in 1:k){
    P.matrix[,j] <- sapply(1:nrow(Z.matrix), 
        function(i) pchisq(Z.matrix[i,j]^2,df=1,ncp=0,lower.tail=F))
}	
colnames(Z.matrix) <- paste("Z",1:k,sep="")
colnames(P.matrix) <- paste("P",1:k,sep="")

## Steps to implementing PLACO+
# Step 1: Obtain the variance and correlation parameters estimates (only once)
VarZ <- var.placo(Z.matrix, P.matrix, p.threshold = 1e-4) # Calculate variance estimates for each trait
CorZ <- cor.pearson(Z.matrix, P.matrix, p.threshold = 1e-4, returnMatrix=F) # Calculate correlation parameter estimate between traits

# Step 2: Apply test of pleiotropy for each variant using PLACO+
out <- sapply(1:p, function(i) placo.plus(Z=Z.matrix[i,], VarZ=VarZ, CorZ=CorZ))
# (returns a 2xp matrix of PLACO+ test statistics and p-values for p variants)

# Check the output for say variant 100
dim(out)
out[,100]$T.placo.plus
out[,100]$p.placo.plus
