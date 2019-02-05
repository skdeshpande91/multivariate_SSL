# Example script
library(SSLASSO)
library(MASS)
library(mSSL)


n <- 400
p <- 500
q <- 25
rho <- 0.9
set.seed(129)
SigmaX <- 0.7^(abs(outer(1:p, 1:p, "-")))
X <- mvrnorm(n, rep(0, times = p), SigmaX)

B.orig <- matrix(sample(c(runif(floor(p*q*0.2), -2,2), rep(0, times = p*q - floor(p*q*0.2)))), nrow = p, ncol = q)
XB <- X %*% B.orig

Sigma.orig <- rho^(abs(outer(1:q, 1:q, "-")))
Omega.orig <- solve(Sigma.orig)
# R's solve() will sometimes not set entries in Omega exactly equal to 0. 
Omega.orig[abs(Omega.orig) < 1e-3] <- 0


L <- 10
lambdas <- list(lambda1 = 1, lambda0 = seq(10, n, length = L))
xis <- list(xi1 = 0.01*n, xi0 = n*seq(0.1, 1, length = L))

set.seed(718)
E <- mvrnorm(n, rep(0, times = q), Sigma.orig)
Y <- XB + E


fit_mSSL_dcpe <- mSSL_dcpe(X, Y)
fit_sep_sslg <- sep_sslg(X,Y)


# Takes substantially longer to run
fit_mSSL_dpe <- mSSL_dpe(X,Y)
