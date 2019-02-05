make.pd <- function (m, cond.num) 
{
    if (!is.matrix(m)) 
        m <- as.matrix(m)
    d <- dim(m)[1]
    if (dim(m)[2] != d) 
        stop("Input matrix is not square!")
    es <- eigen(m)
    esv <- es$values
  
    esv.max <- max(esv)
    if (esv.max<=0)
      stop("The largest eigenvalue is not positive. Try m'=-m.")
    esv.min <- min(esv)
    delta <- (esv.max - cond.num*esv.min)/(cond.num-1)
     
    m2 <- es$vectors %*% diag(esv+delta, d) %*% t(es$vectors)
    return(m2)
}
