###------------- Code for formally testing pleiotropic association of two phenotypes with SNPs using GWAS summary statistics-----------------
#
message("====================================================================")
message("                  PLACO v0.2.0 is loaded")
message("====================================================================")
message("If you use this software, please cite:")
message("")
message("Original PLACO paper -->")
message("Ray and Chatterjee (2020) A powerful method for pleiotropic analysis")
message("    under composite null hypothesis identifies novel shared")
message("    loci between type 2 diabetes and prostate cancer.")
message("    PLoS Genetics 16(12): e1009218")
message("")
message("PLACO+ paper -->")
message("Park and Ray (2025+) A robust pleiotropy method with")
message("    applications to lipid traits and to inflammatory bowel")
message("    disease subtypes with sample overlap. Submitted.")
message("")
message("**See https://github.com/RayDebashree/PLACO for updated citations")
message("--------------------------------------------------------------------")
message("")

############################################
#---------------- Function for normal product based tail probability calculation 
# (Using modified Bessel function of the 2nd kind with order 0)
.pdfx<-function(x) besselK(x=abs(x),nu=0)/pi
# Function for calculating the p-value under the composite null hypothesis using (symmetric) normal product distribution
.p.bessel<-function(z, varz, AbsTol=1e-13){
	p1 <- 2*as.double(integrate(Vectorize(.pdfx), abs(z[1]*z[2]/sqrt(varz[1])), Inf, abs.tol=AbsTol)$value)
	p2 <- 2*as.double(integrate(Vectorize(.pdfx), abs(z[1]*z[2]/sqrt(varz[2])), Inf, abs.tol=AbsTol)$value)
	p0 <- 2*as.double(integrate(Vectorize(.pdfx), abs(z[1]*z[2]), Inf, abs.tol=AbsTol)$value)
	pval.compnull <- p1 + p2 - p0
	return(pval.compnull)
}

#---------------- Function for bivariate normal product with correlation based tail probability calculation 
# (Using modified Bessel function of the 2nd kind with order 0)
.pdfx.cor<-function(x, corz){ 
  # Calculate the modified Bessel function of the second kind for tail probability adjustment
  f1 <- besselK(x=(abs(x)/(1-corz^2)), nu=0)
  f2 <- exp(corz*x/(1-corz^2))/(pi*sqrt(1-corz^2))
  y<-f1*f2
  for(i in 1:length(y)){
    if(f1[i]==0){ y[i]<-0 }else{ y[i] <- f1[i]*f2[i] }
  }
  return(y)}

# Modified function for calculating the p-value under the composite null hypothesis using (asymmetric) bivariate normal product distribution
.p.bessel.cor<-function(z, varz, corz, AbsTol=1e-13){
  p1 <- as.double(integrate(Vectorize(.pdfx.cor), abs(z[1]*z[2]/sqrt(varz[1])), Inf, corz=corz, abs.tol=AbsTol)$value + integrate(Vectorize(.pdfx.cor), -Inf, -abs(z[1]*z[2]/sqrt(varz[1])), corz=corz, abs.tol=AbsTol)$value)
  p2 <- as.double(integrate(Vectorize(.pdfx.cor), abs(z[1]*z[2]/sqrt(varz[2])), Inf, corz=corz, abs.tol=AbsTol)$value + integrate(Vectorize(.pdfx.cor), -Inf, -abs(z[1]*z[2]/sqrt(varz[2])), corz=corz, abs.tol=AbsTol)$value)
  p0 <- as.double(integrate(Vectorize(.pdfx.cor), abs(z[1]*z[2]), Inf, corz=corz, abs.tol=AbsTol)$value + integrate(Vectorize(.pdfx.cor), -Inf, -abs(z[1]*z[2]), corz=corz, abs.tol=AbsTol)$value)
  # Compute final p-value under composite null hypothesis
  pval.compnull <- p1 + p2 - p0
  return(pval.compnull)
}


#---------------- Function for estimating the variances for PLACO
var.placo<-function(Z.matrix, P.matrix, p.threshold=1e-4){
	# Here Z.matrix is the pxk matrix of Z-scores where p is total no. of variants in the dataset, and k is the no. of traits 
	# Similarly, P.matrix is the corresponding pxk matrix of p-values where p is total no. of variants in the dataset, and k is the no. of traits
	# p.threshold determines which variants are deemed to have no association marginally (default: 1e-4)
  # checks
  k<-ncol(Z.matrix)
  if(k!=2) stop("This method is meant for 2 traits only. Columns correspond to traits.")
  ZP<-cbind(Z.matrix,P.matrix)
  ZP<-na.omit(ZP)

  rows.alt<-which(ZP[,3]<p.threshold & ZP[,4]<p.threshold)
  if(length(rows.alt)>0){
    ZP<-ZP[-rows.alt,]
    if(nrow(ZP)==0) stop(paste("No 'null' variant left at p-value threshold",p.threshold))
    if(nrow(ZP)<30) warning(paste("Too few 'null' variants at p-value threshold",p.threshold))
  }
  varz<-diag(var(ZP[,c(1,2)]))
  return(varz)
}

#---------------- Function for estimating correlation matrix of the Z's
cor.pearson<-function(Z.matrix, P.matrix, p.threshold=1e-4, returnMatrix=TRUE){
	# Here Z.matrix is the pxk matrix of Z-scores where p is total no. of variants in the dataset, and k is the no. of traits 
	# Similarly, P.matrix is the corresponding pxk matrix of p-values where p is total no. of variants in the dataset, and k is the no. of traits
	# p.threshold determines which variants are deemed to have no association marginally (default: 1e-4)
  # If returnMatrix=TRUE (default), this function returns the kxk correlation matrix
  # If k=2 and returnMatrix=FALSE, this function returns the correlation parameter estimate ([1,2]th element of correlation matrix)
  # checks
  k <- ncol(Z.matrix)
  if(k!=2) stop("This method is meant for 2 traits only.")
  # estimating correlation
    row.exclude <- which( apply(P.matrix, MARGIN = 1, function(x) any(x < p.threshold)) == TRUE )
    if(length(row.exclude)>0) Z.matrix<-Z.matrix[-row.exclude,]
  R <- cor(Z.matrix)
  if(isTRUE(returnMatrix)){ return(R) }else{ return(R[1,2]) }
}

############################################
#---------------- Function for implementing original PLACO (Ray et al, 2020) - applicable to independent traits only
placo <- function(Z, VarZ, AbsTol=.Machine$double.eps^0.8){
				# Z: vector of Z-scores of size k=2 (i.e., collection of Z-scores of a particular SNP for k=2 traits)
				# VarZ: vector of variances of Z-scores (covariance assumed 0; so need to be independent traits)
				# AbsTol: absolute tolerance (accuracy parameter) for numerical integration.
   # checks				
   k <- length(Z)		
   if(k!=2) stop("This method is meant for 2 traits only.")
   if(length(VarZ)!=k) stop("Provide variance estimates for 2 traits as obtained using var.placo() function.")
   
   # test of pleiotropy: PLACO
   pvalue.b = .p.bessel(z=Z, varz=VarZ, AbsTol=AbsTol)
   return(list(T.placo=prod(Z), p.placo=pvalue.b))
}

#---------------- Function for implementing PLACO+ (Park et al, 2025) - applicable to any two traits 
placo.plus <- function(Z, VarZ, CorZ, AbsTol=.Machine$double.eps^0.8){
      # Z: vector of Z-scores of size k=2 (i.e., collection of Z-scores of a particular SNP for k=2 traits)
      # VarZ: vector of variances of Z-scores 
      # CorZ : estimated correlation between two traits (scalar quantity, not a matrix)
      # AbsTol: absolute tolerance (accuracy parameter) for numerical integration.
  # checks				
  k <- length(Z)		
  if(k!=2) stop("This method is meant for 2 traits only.")
  if(length(VarZ)!=k) stop("Provide variance estimates for 2 traits as obtained using var.placo() function.")
  if(is.matrix(CorZ)) stop("CorZ is estimated correlation parameter (scalar), not a matrix. Use cor.pearson() function with returnMatrix=FALSE to obtain it.")
  if(CorZ < -1 | CorZ > 1) stop("Estimated correlation parameter (CorZ) must be between -1 and 1. Use cor.pearson() function to obtain it.")
  
  # robust, general test of pleiotropy: PLACO+
  pvalue.b = .p.bessel.cor(z=Z, varz=VarZ, corz=CorZ, AbsTol=AbsTol)
  return(list(T.placo.plus=prod(Z), p.placo.plus=pvalue.b))
}



