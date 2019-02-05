capme <- function( x, y, 
                   lambda = NULL,  
                   nlambda = ifelse(is.null(lambda),10,length(lambda)),
                   lambda.max = max(cov(x)*(1-1/nrow(x)))*max(cov(y)*(1-1/nrow(y)))*
                                   sqrt(4*log(ncol(x))*log(ncol(y))/n),
                   lambda.min = ifelse(nrow(x)>ncol(x), 1e-6, 1e-4),
                   logspaced.lambda = TRUE,
                   linsolver.Gamma = c("simplex","primaldual"), 
                   ntau = ifelse(is.null(tau),50,length(tau)),
                   tau = NULL,
                   tau.max = 0.8,
                   tau.min = ifelse(nrow(x)>ncol(x), 1e-6, 1e-4),
                   perturb = TRUE,
                   logspaced.tau = TRUE,
                   linsolver.Omega = c("simplex","primaldual"),
                   pdtol = 1e-3, pdmaxiter = 50)
{                       
  lpfun.Gamma <- match.arg(linsolver.Gamma, c("simplex","primaldual"))
  lpfun.Omega <- match.arg(linsolver.Omega, c("simplex","primaldual"))
  
  n <- nrow(x)
  q <- ncol(x)
  if(nrow(y)!= n){
    stop("The dimensions of x and y do not match! Try different x and y.")
  }
  p <- ncol(y)

  x <- scale(x)
  y <- scale(y,center=TRUE,scale=FALSE)
                       
  if (is.null(lambda)) {
    if (logspaced.lambda) {
      lambda <- 10^(seq(log10(lambda.min), log10(lambda.max), length.out=nlambda))
    } else {
      lambda <- seq(lambda.min, lambda.max, length.out=nlambda)
    }
  }
   
  if (is.null(tau)) {
    if (logspaced.tau) {
      tau <- 10^(seq(log10(tau.min), log10(tau.max), length.out=ntau))
    } else {
      tau <- seq(tau.min, tau.max, length.out=ntau)
    }
  }
                      
  dantzig.obj <- dantzig( x, y, 
                    lambda = lambda,  
                    logspaced = logspaced.lambda,
                    perturb = perturb,
                    linsolver = lpfun.Gamma, 
                    pdtol = 1e-3, pdmaxiter = 50)                 
    
  Gammalist <- dantzig.obj$Gammalist     
  
  Omegalist <- NULL              
  for (k in 1:nlambda){ 
    dan.obj <- dantzig.obj                  
    dan.obj$lambda <- lambda[k]
    dan.obj$Gammalist <- list(dantzig.obj$Gammalist[[k]])  
    meanclime.obj <- meanclime( dan.obj, 
                       tau = tau, 
                       perturb = perturb,
                       logspaced = logspaced.tau,
                       linsolver.Omega = lpfun.Omega,
                       pdtol = 1e-3, pdmaxiter = 50 )
    Omegalist <- c( Omegalist, meanclime.obj$Omegalist )
  }                    
  
  outlist <- list(Gammalist = Gammalist, Omegalist = Omegalist, x = x, y = y,
                  lambda = lambda, tau = tau, perturb = perturb, 
                  lpfun.Gamma = lpfun.Gamma, lpfun.Omega = lpfun.Omega)
  class(outlist) <- c("capme")
  return(outlist)                     
}
