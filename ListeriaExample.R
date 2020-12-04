data(listeria)

# Summaries
summary(listeria)
plotMissing(listeria)
plotMap(listeria)

for (chrom in 1:19) {
  #chrom=10
  #image(listeria$geno[[chrom]]$data)
  #sum_is.na <- function(x) sum(is.na(x))
  #Missing values per marker
  #apply(listeria$geno[[chrom]]$data, 2, sum_is.na)
  #Missing values per individual
  #apply(listeria$geno[[chrom]]$data, 1, sum_is.na)
  
  SNP.data1 <- listeria$geno[[chrom]]$data
  SNP.qtl <- create.QTLcross(SNP.data1, rnorm(nrow(SNP.data1)))
  
  TEST =TRUE
  start.shrink = 1
  while(TEST) {
    map <- get.map(SNP.qtl, shrink=start.shrink)
    delta <- 0.01 + (map$conv1>10)*0.05 + (map$conv1>2)*0.05
    start.shrink <- start.shrink - delta
    TEST <- (map$conv1 + map$conv2) > 0
  }
  #plot(map$indx)
 
  image( Matrix( map$sparse.prec[map$indx, map$indx] ) , main=paste("Chromosome", chrom))  #Should be close to tridiagonal. can be used to assess the fit
  
  cat("Chromosome no.",chrom, "\n")
  check.order(Matrix( map$sparse.prec[map$indx, map$indx] ))
  if (all(map$indx==(1:length(map$indx))) | all(map$indx==seq(length(map$indx),1) )) print("Correct ordering")
      #listeria$geno[[chrom]]$map
      
}
