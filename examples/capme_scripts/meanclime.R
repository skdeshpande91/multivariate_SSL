meanclime <- function(dantzig.obj, 
                  tau=NULL, 
                  ntau=ifelse(is.null(tau),100,length(tau)),
                  tau.max=0.8,
                  tau.min=ifelse(nrow(x)>ncol(x), 1e-3, 1e-2),
                  perturb=FALSE,
                  logspaced=TRUE,
                  linsolver.Omega=c("primaldual", "simplex"),
                  pdtol=1e-3, pdmaxiter=50
                  )
{
  lpfun <- match.arg(linsolver.Omega, c("primaldual", "simplex"))
  lambda <- dantzig.obj$lambda[1]
    
  x <- dantzig.obj$x
  y <- dantzig.obj$y
  Gammalist <- dantzig.obj$Gammalist
  if( length(Gammalist)>1 ){
    warning("Too many Gammas in the Gammalist of dantzig.obj. Only the first one is used.")
  }
  Gamma <- Gammalist[[1]]
  re <- y - x%*%Gamma
  n <- nrow(x)
  q <- ncol(x)
  p <- ncol(y)
  
    
  if (is.null(tau)) {
    if (logspaced) {
      tau <- 10^(seq(log10(tau.min), log10(tau.max), length.out=ntau))
    } else {
      tau <- seq(tau.min, tau.max, length.out=ntau)
    }
  }
  
  Sigma <- var(re)*(1-1/n)  
  
  ## Set to perturbed Sigma to have conditional number p
  eigvals <- eigen(Sigma, only.values=T)$values
  if (is.logical(perturb)) {
      if (perturb) { 
          perturb <- max(max(eigvals) - p*min(eigvals), 0)/(p-1)
      } else {
          perturb <- 0
      }
  }
  
  Sigma <- Sigma+1.5*diag(p)*perturb
  emat <- diag(p)
  
  Omegalist <- vector("list", ntau)
  if (lpfun == "simplex") {
    for (j1 in 1:ntau) {
      Omega <- matrix(0, nrow=p, ncol=p)
      tau.j1 <- tau[j1]
      for (j in 1:p) {
        beta <- linprogS2(Sigma, emat[,j], tau.j1)
        Omega[,j] <- beta
      }
      Omegalist[[j1]] <- Omega*(abs(Omega)<= abs(t(Omega)))+ t(Omega)*(abs(Omega)> abs(t(Omega)))
    }
  }
  
  if (lpfun == "primaldual") {
    Omega0 <- solve(Sigma)
    
    for (j1 in 1:ntau) {
      Omega <- matrix(0, nrow=p, ncol=p)
      tau.j1 <- tau[j1]
      for (j in 1:p) {
        #cat("j,j1", j,j1,"\n")
        beta <- linprogPD2(Omega0[,j], Sigma, emat[,j], tau.j1, pdtol, pdmaxiter)
        Omega[,j] <- beta
      }
      Omegalist[[j1]] <- Omega*(abs(Omega)<= abs(t(Omega)))+ t(Omega)*(abs(Omega)> abs(t(Omega)))
    }
  }


  outlist <- list(Gamma=Gamma, Omegalist=Omegalist, x = x, y = y, lambda = lambda, 
                  tau = tau, perturb=perturb, 
                  lpfun.Gamma=dantzig.obj$lpfun, lpfun.Omega=lpfun)
  class(outlist) <- c("meanclime")
  return(outlist)
}
