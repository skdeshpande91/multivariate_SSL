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

B_orig <- matrix(sample(c(runif(floor(p*q*0.2), -2,2), rep(0, times = p*q - floor(p*q*0.2)))), 
                 nrow = p, ncol = q)
XB <- X %*% B_orig

Sigma_orig <- rho^(abs(outer(1:q, 1:q, "-")))
Omega_orig <- solve(Sigma_orig)
# R's solve() will sometimes not set entries in Omega exactly equal to 0. 
Omega_orig[abs(Omega_orig) < 1e-3] <- 0


set.seed(724)
E <- mvrnorm(n, rep(0, times = q), Sigma_orig)
Y <- XB + E

# May need to run it with a finer grid. It seems like DPE over-selected Beta and has diagonal Omega.
fit_mSSL_dcpe <- mSSL_dcpe(X,Y)

fit_mSSL_dpe <- mSSL_dpe(X,Y)

save(fit_mSSL_dcpe, fit_mSSL_dpe, X, Y,E, B_orig, Omega_orig, file = "feb7.RData")