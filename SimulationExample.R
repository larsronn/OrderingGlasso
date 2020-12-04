############### SIMULATION FUNCTION #################
simulate.cross <- function(LD, n.ind, m.markers, n.chr = 1, err = 1e-6, remove.identical = FALSE) {
  n <- n.ind
  m <- m.markers
  SNP1 <- SNP2 <- matrix(0, n, m*n.chr)
  for (i.chr in 1:n.chr) {
    for (j in 1:m) {
      if (j==1) {
        vec1 <- rbinom(n, 1, 0.5)
        vec2 <- rbinom(n, 1, 0.5)
      }
      if (j>1) vec1 <- abs(vec1-rbinom(n, 1, LD))
      if (j>1) vec2 <- abs(vec2-rbinom(n, 1, LD))
      SNP1[ ,(j+(i.chr-1)*m)] <- vec1
      SNP2[ ,(j+(i.chr-1)*m)] <- vec2
    }
  }
  #Add some noise specified by err 
  SNP1 <- abs(SNP1 - matrix(rbinom(n*ncol(SNP1), 1, err), n, ncol(SNP1)))
  SNP2 <- abs(SNP2 - matrix(rbinom(n*ncol(SNP2), 1, err), n, ncol(SNP2)))
  SNP <- SNP1 + SNP2 #Diploid
  if (remove.identical) {
    to.remove = NULL
    for (i in 2:m.markers) {
      if (abs(cor(SNP[,i],SNP[,i-1]))==1) to.remove <- c(to.remove, i)
    }
    if (!is.null(to.remove)) SNP <- SNP[,-to.remove]
  }
  #Output matrix with elements having values 1, 2, or 3
  return( SNP + 1 )
}

############### FUNCTION MAKING SIMULATED DATA INTO R/qtl OBJECT #############
create.QTLcross <- function(SNP.data, phen) {
	N <- nrow(SNP.data)
	m <- ncol(SNP.data)
	if (length(phen)!=N) stop("phenotype vector of wrong length")
	tmp1 <- NULL
	X1 <- as.integer(SNP.data)
	dim(X1) <- c(N,m)
	colnames(X1) <- as.character(1:m)
	tmp1$geno[[1]]$data <- X1
	map1 <- as.numeric(1:m)
	names(map1) <- as.character(1:m)
	tmp1$geno[[1]]$map <- map1
	names(tmp1$geno)="1"
	class(tmp1$geno[[1]])="A"
	tmp1$pheno <- data.frame(y=phen)
	attr(tmp1, "class") <- c("f2", "cross")
	return(tmp1)
}


########### MAIN CODE #################
run.example=TRUE
if (run.example) {
  
  set.seed(1234)
  n=50   #Number of simulated individuals
  m=100   #Number of markers per chromosome
  SNP <- simulate.cross(LD = 0.1, n.ind = n, m.markers = m, n.chr = 1 )
  reord <- sample(1:ncol(SNP), ncol(SNP) )
  SNP <- SNP[,reord] #Shuffle order
  SNP.qtl <- create.QTLcross(SNP, phen = rnorm(n))
  TEST =TRUE
  start.shrink = 1
  while(TEST) {
    map <- get.map(SNP.qtl, shrink=start.shrink)
    delta <- 0.01 + (map$conv1>10)*0.05 + (map$conv1>2)*0.05
    start.shrink <- start.shrink - delta
    TEST <- (map$conv1 + map$conv2) > 0
  }
  
  image( Matrix( map$sparse.prec[map$indx, map$indx] ) ) #Should be close to tridiagonal. can be used to assess the fit
  check.order( Matrix( map$sparse.prec[map$indx, map$indx] ) )
  plot( reord[map$indx] , xlab="True Order", ylab="Estimated Order") #Compare true and estimated ordering
  
}

