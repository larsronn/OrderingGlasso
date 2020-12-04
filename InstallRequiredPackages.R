if (!("INLA" %in% installed.packages())) {
  install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
}
if (!("glasso" %in% installed.packages())) {
  install.packages("glasso")
}
