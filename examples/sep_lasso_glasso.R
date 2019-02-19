######################
# Separate LASSO followed by GLASSO
######################

library(glmnet)
library(glasso)
library(cvTools)
################
# Need to write our own cv.glasso
# Essentially use this: https://github.com/czarrar/causality_primer/blob/master/journal/2014-07-12_cross_validation_graphical_lasso.Rmd

log.likelihood <- function(Omega, S){
  q <- ncol(Omega)
  log.det <- determinant(Omega, logarithm = TRUE)$modulus
  log.like <- log.det - sum(diag(S %*% Omega))
  return(log.like)
}


cv.glasso <- function(R, n){
  S.full <- t(R) %*% R/n
  full.fit <- glassopath(S.full, trace = 0)
  xi.list <- full.fit$rholist
  
  folds <- cvFolds(n,10, type = "consecutive")
  log.like <- sapply(1:10, function(ki){
    train.index <- which(folds$which != ki)
    test.index <- which(folds$which == ki)
    n.train <- length(train.index)
    n.test <- length(test.index)
    S.train <- t(R[train.index,]) %*%R[train.index,]/n.train
    S.test <- t(R[test.index,]) %*% R[test.index,]/n.test
    GLP <- glassopath(S.train, xi.list, trace = 0)
    tmp <- apply(GLP$wi, 3, function(Omega.hat) log.likelihood(Omega.hat, S.test))
    tmp
  })
  # returns a 10 x 10 matrix, columns correspond to number of folds, rows for a specific value of xi
  max.ind <- which.max(rowMeans(log.like))
  return(full.fit$wi[,,max.ind])
 

}

sep.lasso.glasso <- function(Y, X){

  n <- nrow(Y)
  q <- ncol(Y)
  p <- ncol(X)
  
  # No intercept
  B.hat <- matrix(nrow = p, ncol = q)
  for(k in 1:q){
    #tmp.fit <- cv.glmnet(X, Y[,k], family = "gaussian", intercept = FALSE)
    tmp.fit <- cv.glmnet(X, Y[,k], family = "gaussian", intercept = TRUE)
    tmp.bhat <- as.vector(coef(tmp.fit, method = "lambda.1se"))[-1] # drop the intercept
    B.hat[,k] <- tmp.bhat
 
    
  }
  R <- Y - X %*% B.hat
  Omega.hat <- cv.glasso(R,n)
  # Do the glasso now
  return(list(B = B.hat, Omega = Omega.hat))
}