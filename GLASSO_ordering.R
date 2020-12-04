library(qtl) # R/qtl package used for to computing recombination fractions
library(Matrix)
library(glasso) #Graphical lasso function
#INLA ackage used to access the ordering function
#install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
library(INLA) 

############ FUNCTION USING GRAPHICAL LASSO PLUS RE-ORDERING FUNCTION TO FIND TRUE ORDER #############
get.map <- function(qtl.obj, shrink=0.9, verb=FALSE) {
  test <- est.rf(qtl.obj)
  RF <- t(test$rf)
  RF[lower.tri(RF)] <- t(RF)[lower.tri(RF)]
  diag(RF) <- 0
  KORR <- 1- 2*RF #Transform recombination fractions to korrelations
  gCM <- glasso(KORR, rho=shrink)
  test1 <- inla.qreordering(gCM$wi, reordering="band")
  indx <- test1$ireordering
  conv <- sum(colSums(gCM$wi[indx,indx]) - diag(gCM$wi[indx,indx]) ==0)
  if (verb) {
    cat("Number of uncorrelated markers", conv, "\n")
    cat("Number of independent clusters", sum( colSums(tril(gCM$wi[indx, indx], -1)!=0) == 0) - 1 , "\n" )
  }
  #cat( colSums(tril(gCM$wi[indx,indx], -1)!=0), "\n" )
  list(indx = indx, sparse.prec = gCM$wi, conv1 = conv, conv2 = sum( colSums(tril(gCM$wi[indx, indx], -1)!=0) == 0) - 1)
}

check.order <- function(M) {
  for (i in 1:ncol(M)) {
    x <- M[i:nrow(M),i]
    y <- x[x!=0]
    if (!all(order(abs(y)) == seq(length(y), 1))) print("Warning: Matrix not in perfect order")
  }
}


