cv.capme <- function(fold = 5, loss = c("likelihood", "tracel2"), x, y, 
                     lambda = NULL,  
                     nlambda = ifelse(is.null(lambda),10,length(lambda)),
                     lambda.max = max(cov(x)*(1-1/nrow(x)))*max(cov(y)*(1-1/nrow(y)))*
                                   sqrt(4*log(ncol(x))*log(ncol(y))/n),
                     lambda.min=ifelse(nrow(x)>ncol(x), 1e-6, 1e-4),
                     logspaced.lambda = TRUE,
                     linsolver.Gamma = c("simplex","primaldual"), 
                     ntau = ifelse(is.null(tau),50,length(tau)),
                     tau = NULL,
                     tau.max = 0.8,
                     tau.min = ifelse(nrow(x)>ncol(x), 1e-6, 1e-4),
                     perturb = TRUE,
                     logspaced.tau = TRUE,
                     linsolver.Omega = c("simplex","primaldual"))
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
      lambda <- 10^(seq(log10(lambda.min), log10(lambda.max), length.out = nlambda))
    } else {
      lambda <- seq(lambda.min, lambda.max, length.out=nlambda)
    }
  }
   
  if (is.null(tau)) {
    if (logspaced.tau) {
      tau <- 10^(seq(log10(tau.min), log10(tau.max), length.out = ntau))
    } else {
      tau <- seq(tau.min, tau.max, length.out = ntau)
    }
  }
        
  part.list <- cv.part(n, fold)
  
  lossname <- match.arg(loss, c("likelihood", "tracel2"))
  lossfun <- match.fun(lossname )
  
  ntau <- length(tau)
  nlambda <- length(lambda)
  
  loss.re <- array(0, c(nlambda,ntau,fold))
  
  for (k in 1:fold){
    x.train <- x[part.list$trainMat[,k],]
    y.train <- y[part.list$trainMat[,k],]
    capme.obj <- capme( x = x.train, y = y.train, lambda = lambda, tau = tau, 
                      linsolver.Gamma = lpfun.Gamma, linsolver.Omega = lpfun.Omega,
                      perturb = perturb)
    Gammalist <- capme.obj$Gammalist
    Omegalist <- capme.obj$Omegalist
    x.test <- x[part.list$testMat[,k],]
    y.test <- y[part.list$testMat[,k],]
    ntest <- nrow(x.test)
    for (i in 1:nlambda){
      for (j in 1:ntau){
        idx <- (i-1)*nlambda + j
        Gamma <- Gammalist[[i]]
        Omega <- Omegalist[[idx]]
        re <- y.test - x.test %*% Gamma
        S <- (1-1/ntest)*var(re)
        loss.re[i,j,k] <- lossfun(S,Omega)
      }
    } 
  }
  
  loss.mean <- apply(loss.re, 1:2, mean)
  loss.sd <- apply(loss.re, 1:2, sd)
  
  ind <- which.min(loss.mean)
  tauind <- ceiling(ind/nlambda)
  lambdaind <- ind - (tauind-1)*nlambda

  tauopt <- tau[tauind]
  lambdaopt <- lambda[lambdaind]
  
 
  outlist <- list(tauopt = tauopt, tauind = tauind, lambdaopt = lambdaopt,
                  lambdaind = lambdaind, loss = lossname, tau = tau, lambda = lambda,
                  loss.mean = loss.mean, loss.sd = loss.sd)
  class(outlist) <- c("cv.capme")
  return(outlist) 
}
