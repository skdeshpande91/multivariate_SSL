Mpower <- function(Sig,q)
{
    a <- svd(Sig);
    d <- a$d^(q);
    m <- a$u%*%diag(d)%*%t(a$v);
    m   
}
