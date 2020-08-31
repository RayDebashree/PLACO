###------------- Code for formally testing pleiotropic association of two phenotypes with SNPs using GWAS summary statistics-----------------
#
message("============================================================")
message("                  PLACO v0.1.1 is loaded")
message("============================================================")
message("If you use this software, please cite:")
message("Ray et al.(2020) A powerful method for pleiotropic analysis")
message("    under composite null hypothesis identifies novel shared")
message("    loci between type 2 diabetes and prostate cancer.")
message("    BioRxiv https://doi.org/10.1101/2020.04.11.037630")
message("------------------------------------------------------------")
message("")

############################################
#---------------- Function for normal product based tail probability calculation 
# (Using modified Bessel function of the 2nd kind with order 0)
.pdfx<-function(x) besselK(x=abs(x),nu=0)/pi
.p.bessel<-function(z, varz, AbsTol=1e-13){
	p1<-2*as.double(integrate(Vectorize(.pdfx), abs(z[1]*z[2]/sqrt(varz[1])),Inf, abs.tol=AbsTol)$value)
	p2<-2*as.double(integrate(Vectorize(.pdfx), abs(z[1]*z[2]/sqrt(varz[2])),Inf, abs.tol=AbsTol)$value)
	p0<-2*as.double(integrate(Vectorize(.pdfx), abs(z[1]*z[2]),Inf, abs.tol=AbsTol)$value)
	pval.compnull<-p1+p2-p0
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
cor.pearson<-function(Z.matrix, P.matrix, p.threshold=1e-4){
	# Here Z.matrix is the pxk matrix of Z-scores where p is total no. of variants in the dataset, and k is the no. of traits 
	# Similarly, P.matrix is the corresponding pxk matrix of p-values where p is total no. of variants in the dataset, and k is the no. of traits
	# p.threshold determines which variants are deemed to have no association marginally (default: 1e-4)
  # checks
  k<-ncol(Z.matrix)
  if(k!=2) stop("This method is meant for 2 traits only.")
  # estimating correlation
    row.exclude<-which( apply(P.matrix, MARGIN = 1, function(x) any(x < p.threshold)) == TRUE )
    if(length(row.exclude)>0) Z.matrix<-Z.matrix[-row.exclude,]
  R<-cor(Z.matrix)
  return(R)
}

############################################
placo<-function(Z, VarZ, AbsTol=.Machine$double.eps^0.8){
				# Z: vector of Z-scores of size k=2 (i.e., collection of Z-scores of a particular SNP for k=2 traits)
				# VarZ: vector of variances of Z-scores (covariance assumed 0; so need to be independent traits)
				# AbsTol: absolute tolerance (accuracy paramater) for numerical integration.
   # checks				
   k<-length(Z)		
   if(k!=2) stop("This method is meant for 2 traits only.")
   if(length(VarZ)!=k) stop("Provide variance estimates for 2 traits as obtained using var.placo() function.")
   
   # test of pleiotropy: PLACO
   pvalue.b=.p.bessel(z=Z, varz=VarZ, AbsTol=AbsTol)
   return(list(T.placo=prod(Z), p.placo=pvalue.b))
}
