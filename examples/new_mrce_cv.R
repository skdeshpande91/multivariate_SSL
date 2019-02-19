# new MRCE code

# Code taken from the package itself. This is the internal stuff used inside of mrce()

rblasso <-function(s, m, om, nlam, n, B0=NULL, soft=NULL, objective=0, tol=1e-5, maxit=500, quiet=TRUE)
{
  p=dim(s)[1]
  q=dim(om)[1]
  if(is.null(B0))
    B0=as.double(rep(0, p*q))
  if(is.null(soft))
   soft=as.double(m)  
 
  objective=as.double(objective)
  soft=as.double(soft)
  s = as.double(s)
  m = as.double(m)
  om = as.double(om)
  nlam=as.double(nlam)
  tol=as.double(tol)
  totalit=0
  mode(n) = "integer"
  mode(p) = "integer"
  mode(q) = "integer"
  mode(maxit) = "integer"
  mode(totalit)="integer"
  dotCoutput=.C("blasso", B=B0, S=s, M=m, Om=om, soft=soft, pin=p,
            qin=q, nin=n, lam=nlam, tol=tol, maxit=maxit, totalit=totalit, objective=objective)
	
  if(!quiet)
  {  
    cat("Total iterations for solving for B was", dotCoutput$totalit, "\n")
  }
  B = matrix(dotCoutput$B, nrow=p, ncol=q)
  return(list(B=B, obj=dotCoutput$objective))
}

compute.mrce=function(X,Y, lam1, lam2, tol.out, tol.in, maxit.out, maxit.in, silent,
                      cov.tol, cov.maxit,informed=NULL, eps)
{
  n=nrow(X)
  p=ncol(X)
  q=ncol(Y)
  if(is.null(informed))
  {
    if(!is.matrix(lam2))
      nlam=matrix(n*lam2, nrow=p, ncol=q) else nlam=n*lam2
    
    mx=apply(X, 2, mean)
    my=apply(Y, 2, mean)
    X=scale(X, center=mx, scale=FALSE)
    Y=scale(Y, center=my, scale=FALSE)
    yty=crossprod(Y)
    xty=crossprod(X,Y)
    xtx=crossprod(X)
    old.B=matrix(0, nrow=p, ncol=q)
    tolmult=sum(diag(yty)/n)
    tout=tol.out*tolmult
    residual.cov = yty/n
    sigma=diag(diag(residual.cov))
    om=diag(1/diag(residual.cov))
    omoff=om
    diag(omoff)=0
    old.obj=sum(residual.cov*om)-determinant(om, logarithm=TRUE)$mod[1] + lam1*sum(abs(omoff))+2*sum(lam2*abs(old.B))
  } else
  {
    nlam=matrix(n*lam2, nrow=p, ncol=q)
    mx=informed$mx
    my=informed$my
    xtx=informed$xtx
    xty=informed$xty
    yty=informed$yty
    old.B=informed$Bhat
    om=informed$omega
    sigma=informed$sigma
    tolmult=sum(diag(yty)/n)
    tout=tol.out*tolmult
    residual.cov = crossprod(Y-X%*%old.B)/n
    omoff=om
    diag(omoff)=0
    old.obj=sum(residual.cov*om)-determinant(om, logarithm=TRUE)$mod[1] + lam1*sum(abs(omoff))+2*sum(lam2*abs(old.B))
  }
  Moff=lam1*(1-diag(q))
  k=0  
  iterating=TRUE
  while(iterating)
  {
    k=k+1  
    if(min(diag(residual.cov)) < eps)
    {
      cat("A perfect fit occured (lam2 may be too small). Terminated early.\n")
      break
    }   
    if(k == 1)
    {
      cov.out=NULL
      cov.out$X=om
      cov.out$W=sigma
    }
    cov.out=QUIC(S=residual.cov, rho=Moff,tol=cov.tol, msg=0, maxIter=cov.maxit, X.init=cov.out$X, W.init=cov.out$W)
    if(!is.matrix(cov.out$X))
      cov.out$X=matrix(cov.out$X, nrow=q, ncol=q)
    if(!is.matrix(cov.out$W))
      cov.out$W=matrix(cov.out$W, nrow=q, ncol=q)  
    om=cov.out$X   
    tolinmult=sum(yty*om)/n 
    ## added
    omoff=om
    diag(omoff)=0
    obj.after.omega=sum(residual.cov*om)-determinant(om, logarithm=TRUE)$mod[1] + lam1*sum(abs(omoff))+2*sum(lam2*abs(old.B))
    if(!silent) cat("k =", k, "obj. fn. val. after Omega update is", obj.after.omega, "\n")  
	xtyom=xty%*%om  
    soft=xtyom - xtx%*%old.B%*%om + old.B*tcrossprod(diag(xtx), diag(om)) 
    outlasso=rblasso(s=xtx, m=xtyom, om=om, nlam=nlam, n=n,B0=old.B, soft=soft, objective=obj.after.omega, tol=(tol.in*tolinmult), maxit=maxit.in, quiet=silent)		
    old.B=outlasso$B
    residual.cov = crossprod(Y-X%*%old.B)/n
    new.obj=outlasso$obj
    bdist = old.obj-new.obj
    iterating= (bdist > tout) & (k <= maxit.out)
    old.obj=new.obj
    if(!silent) cat("k =", k, "obj. fn. val. after B update is", new.obj, "\n")
  }
  if(!silent) cat("Total outer iterations for MRCE : ", k, "\n")
  muhat=as.numeric(my - crossprod(old.B, mx))
  return(list(Bhat=old.B, muhat=muhat, omega=om, sigma=cov.out$W, mx=mx, my=my))
}


new.mrce.cv <- function(X, Y, lam.vec.1, lam.vec.2, kfolds = 5, tol.out = 1e-8, tol.in = 1e-8, maxit.out = 1e3, maxit.in = 1e3, silent = TRUE, cov.tol = 1e-4,cov.maxit=1e3, eps = 1e-5){
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  
  #err.mat <- matrix(0, nrow = length(lam.vec.1), ncol = length(lam.vec.2))
  
  err.mat.array <- array(NA, dim = c(length(lam.vec.1), length(lam.vec.2), kfolds))
  
  ind <- sample(n)
  #first <- NULL
  for(k in 1:kfolds){
    foldind <- ind[(1 + floor((k-1)*n/kfolds)):floor(k*n/kfolds)]
    X.tr <- X[-foldind,]
    Y.tr <- Y[-foldind,]
    X.va <- X[foldind,,drop=FALSE]
    Y.va <- Y[foldind,,drop=FALSE]
    
    mtrx <- apply(X.tr, 2, mean)
    X.tr <- scale(X.tr, scale = FALSE, center = mtrx)
    X.va <- scale(X.va, scale = FALSE, center = mtrx)
    mtry <- apply(Y.tr, 2, mean)
    Y.tr <- scale(Y.tr, scale = FALSE, center = mtry)
    Y.va <- scale(Y.va, scale = FALSE, center = mtry)
    
    n.tr <- nrow(Y.tr)
    informed <- NULL
    informed$mx <- mtrx
    informed$my <- mtry
    informed$yty <- crossprod(Y.tr)
    informed$xtx <- crossprod(X.tr)
    informed$xty <- crossprod(X.tr, Y.tr)
    informed$Bhat <- matrix(0, nrow=p, ncol=q)
    informed$omega <- diag(n.tr/diag(informed$yty))
    informed$sigma <- diag(diag(informed$yty)/n.tr)
    
	  for(i in 1:length(lam.vec.1)){
      if(i > 1) informed=first
      for(j in 1:length(lam.vec.2)){
        tmp.out <- tryCatch({compute.mrce(X=X.tr,Y=Y.tr, lam1=lam.vec.1[i], lam2=lam.vec.2[j], 
                           tol.out=tol.out, tol.in=tol.in, maxit.out=maxit.out, 
                           maxit.in=maxit.in, silent=silent,
                           cov.tol=cov.tol, cov.maxit=cov.maxit, 
                           informed=informed, eps=eps)},
                           error = function(e) NULL)
        if(!is.null(tmp.out)){
          err.mat.array[i,j,k] <- mean( (Y.va - X.va %*% tmp.out$Bhat)^2)
          informed$Bhat <- tmp.out$Bhat
          informed$omega <- tmp.out$omega
          informed$sigma <- tmp.out$sigma
          if(j == 1) first <- informed
        } else { # re-set informed
          informed <- NULL
          informed$mx <- mtrx
          informed$my <- mtry
          informed$yty <- crossprod(Y.tr)
          informed$xtx <- crossprod(X.tr)
          informed$xty <- crossprod(X.tr, Y.tr)
          informed$Bhat <- matrix(0, nrow=p, ncol=q)
          informed$omega <- diag(n.tr/diag(informed$yty))
          informed$sigma <- diag(diag(informed$yty)/n.tr)
        }
        #if(!is.null(tmp.out)){
        #  err.mat[i,j] <- err.mat[i,j] + mean( (Y.va - X.va %*% tmp.out$Bhat)^2)
        #  err.mat.array[i,j,k] <- mean( (Y.va - X.va %*% tmp.out$Bhat)^2)
        #  informed$Bhat <- tmp.out$Bhat
        #  informed$omega <- tmp.out$omega
        #  informed$sigma <- tmp.out$sigma
        #  if(j == 1) first <- informed
        #}                            
        #tmp.out=compute.mrce(X=X.tr,Y=Y.tr, lam1=lam.vec.1[i], lam2=lam.vec.2[j], 
        #                   tol.out=tol.out, tol.in=tol.in, maxit.out=maxit.out, 
        #                   maxit.in=maxit.in, silent=silent,
        #                   cov.tol=cov.tol, cov.maxit=cov.maxit, 
        #                   informed=informed, eps=eps) # need to wrap this in a try-catch statement
        #err.mat[i,j]=err.mat[i,j]+mean((Y.va-X.va%*%tmp.out$Bhat)^2)
        #err.mat.array[i,j,k] <- mean( (Y.va - X.va %*% tmp.out$Bhat)^2)
        #informed$Bhat=tmp.out$Bhat
        #informed$omega=tmp.out$omega
        #informed$sigma=tmp.out$sigma  
        #if(j==1) first=informed        
      }
    }
  }
  ## assemble the err.mat matrix
  err.mat <- matrix(NA, nrow = length(lam.vec.1), ncol = length(lam.vec.2))
  for(i in 1:length(lam.vec.1)){
    for(j in 1:length(lam.vec.2)){
      tmp <- err.mat.array[i,j,]
      if(!all(is.na(tmp))){
        err.mat[i,j] <- sum(tmp, na.rm = TRUE)
      }
    }
  }
  if(all(is.na(err.mat))){
    return(NULL)
  } else{
    # find the (i,j) for the minimum of err.mat
    tmp <- which.min(err.mat) %% (dim(err.mat)[1])
    tmp <- (tmp != 0)*tmp + (tmp == 0) * (dim(err.mat)[1])
    best.i <- tmp
    best.j <- which.min(err.mat[tmp,])
    best.lam1 <- lam.vec.1[best.i]
    best.lam2 <- lam.vec.2[best.j]
    out <- compute.mrce(X=X,Y=Y, lam1=best.lam1, lam2=best.lam2, 
                        tol.out=tol.out, tol.in=tol.in, maxit.out=maxit.out, 
                        maxit.in=maxit.in, silent=silent,
                        cov.tol=cov.tol, cov.maxit=cov.maxit, 
                        informed=NULL, eps=eps)
  return(list(Bhat=out$Bhat, muhat=out$muhat, omega=out$omega, mx=out$mx, my=out$my, 
         best.lam1=best.lam1, best.lam2=best.lam2, cv.err=err.mat, cv.err.mat.array = err.mat.array))
  }

  



}