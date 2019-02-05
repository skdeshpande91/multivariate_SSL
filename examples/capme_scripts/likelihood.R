likelihood <- function(Sigma, Omega) {
  tmp <- (sum(diag(Sigma%*%Omega))  - log(det(Omega)) - dim(Omega)[1])
  if(is.finite(tmp)) {
    return(tmp)
  } else {
    tmp <- Inf
  }
    return(tmp)
}
