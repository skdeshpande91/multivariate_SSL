dantzig <- function ( x, y, lambda = NULL, nlambda =
                       ifelse(is.null(lambda),100,length(lambda)),
                       lambda.max = max(cov(x)*(1-1/nrow(x)))*max(cov(y)*
                       (1-1/nrow(y)))*sqrt(4*log(ncol(x))*log(ncol(y))/n),
                       lambda.min = ifelse(nrow(x)>ncol(x), 1e-3, 1e-2),
                       logspaced = TRUE, perturb = TRUE,
                       linsolver = c("simplex","primaldual"),
                       pdtol = 1e-3, pdmaxiter = 50 )
{
   lpfun <-match.arg(linsolver, c("simplex","primaldual"))
  
   n <- nrow(x)
   q <- ncol(x)
   if(nrow(y)!= n){
     stop("The dimensions of x and y do not match! Try different x and y.")
   }
   p <- ncol(y)
   
   xn <- var(x)*(1-1/n)
   yn <- var(y)*(1-1/n)
   xyn <- cov(x,y)*(1-1/n)
    
   if (is.null(lambda)) {
      if (logspaced) {
        lambda <- 10^(seq(log10(lambda.min), log10(lambda.max), length.out=nlambda))
      } else {
        lambda <- seq(lambda.min, lambda.max, length.out=nlambda)
      }
   }
   
   Gammalist <- vector("list", nlambda)
   
   if (lpfun == "simplex") {
     for (j1 in 1:nlambda) {
       lam <- lambda[j1]
       Gamma<- matrix(0,nrow=q, ncol=p)
       for (j in 1:p) {
         Gamma[,j] <- linprogS2(xn, xyn[,j], lam)
       }
       Gammalist[[j1]] <- Gamma 
     }
   }
   
   if (lpfun == "primaldual") {
    
     if(n>q){
       Gamma0 <- solve(t(x)%*%x)%*%t(x)%*%y    
     }else{
       if(perturb == TRUE){
           tmp <- x%*%t(x) 
           eigvals <- eigen(tmp, only.values=T)$values
           perturb <- max(max(eigvals) - p*min(eigvals), 0)/(p-1)
       }else{
         tmp <- x%*%t(x)                
         tmp <- tmp + diag(n)*perturb
       }
       Gamma0 <- t(x)%*%solve(tmp)%*%y
     }
     for (j1 in 1:nlambda){
       lam <- lambda[j1]
       Gamma<- matrix(0,nrow=q, ncol=p)
       for (j in 1:p) {
         #cat("Gamma,j,j1",j,j1,"\n")
         Gamma[,j] <- linprogPD2(Gamma0[,j], xn, xyn[,j], lam, pdtol, pdmaxiter)
       }
       Gammalist[[j1]] <- Gamma   
     }
     
   }  
    
  outlist <- list(Gammalist = Gammalist, 
                  x = x, y = y, lambda = lambda, lpfun = lpfun)
  class(outlist) <- c("dantzig")
  return(outlist)

}
