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

fit_dpe <- mSSL_dpe(X,Y)
fit_dcpe <- mSSL_dcpe(X,Y)

#save(fit__dcpe, fit_dpe, X, Y,E, B_orig, Omega_orig, file = "feb8.RData")


# Make a histogram

png("B_histogram.png", width = 8, height = 4, units = "in", res = 600)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0), cex.main = 0.9, cex.lab = 0.9, cex.axis = 0.9, mfrow = c(1,2))
# mSSL-DPE
hist(B_orig[fit_dpe$B != 0 & B_orig != 0], breaks = seq(-2, 2, length = 50), 
     ylim = c(0, 85), col  = rgb(0,0,1,1/3), freq = TRUE, main = "mSSL-DPE", xlab = expression(beta[j*k]))
hist(B_orig[fit_dpe$B == 0 & B_orig != 0], breaks = seq(-2,2,length = 50),
     ylim = c(0, 85), col = rgb(1,0,0,1/3), freq = TRUE, add = TRUE)
legend(x = -2, y = 80, pch = 15, col = c(rgb(0,0,1,1/3), rgb(1,0,0,1/3)), legend = c("True Positive", "False Negative"),
       bty = "n", horiz = TRUE, cex = 0.9)
box()
hist(B_orig[fit_dcpe$B != 0 & B_orig != 0], breaks = seq(-2, 2, length = 50), 
     ylim = c(0, 85), col  = rgb(0,0,1,1/3), freq = TRUE, main = "mSSL-DCPE", xlab = expression(beta[j*k]))
hist(B_orig[fit_dcpe$B == 0 & B_orig != 0], breaks = seq(-2,2,length = 50),
     ylim = c(0, 85), col = rgb(1,0,0,1/3), freq = TRUE, add = TRUE)
legend(x = -2, y = 80, pch = 15, col = c(rgb(0,0,1,1/3), rgb(1,0,0,1/3)), legend = c("True Positive", "False Negative"),
       bty = "n", horiz = TRUE, cex = 0.9)
box()
dev.off()


png("B_histogram_bw.png", width = 8, height = 4, units = "in", res = 600)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0), cex.main = 0.9, cex.lab = 0.9, cex.axis = 0.9, mfrow = c(1,2))
# mSSL-DPE
hist(B_orig[fit_dpe$B != 0 & B_orig != 0], breaks = seq(-2, 2, length = 50), 
     ylim = c(0, 85), col  = rgb(204/255,204/255,204/255,1/2), freq = TRUE, main = "mSSL-DPE", xlab = expression(beta[j*k]))
hist(B_orig[fit_dpe$B == 0 & B_orig != 0], breaks = seq(-2,2,length = 50),
     ylim = c(0, 85), col = rgb(55/255,51/255,51/255,1/2), freq = TRUE, add = TRUE)
legend(x = -2, y = 80, pch = c(15,15), col = c(rgb(204/255,204/255,204/255,1), rgb(51/255,51/255,51/255,1/2)), legend = c("True Positive", "False Negative"),
       bty = "n", horiz = TRUE, cex = 0.9)
box()
hist(B_orig[fit_dcpe$B != 0 & B_orig != 0], breaks = seq(-2, 2, length = 50), 
     ylim = c(0, 85), col  = rgb(204/255,204/255,204/255,1/2), freq = TRUE, main = "mSSL-DCPE", xlab = expression(beta[j*k]))
hist(B_orig[fit_dcpe$B == 0 & B_orig != 0], breaks = seq(-2,2,length = 50),
     ylim = c(0, 85), col = rgb(55/255,51/255,51/255,1/2), freq = TRUE, add = TRUE)
legend(x = -2, y = 80, pch = c(15,15), col = c(rgb(204/255,204/255,204/255,1), rgb(51/255,51/255,51/255,1/2)), legend = c("True Positive", "False Negative"),
       bty = "n", horiz = TRUE, cex = 0.9)
box()
dev.off()



# Make a support stability plot
png("support_trajectory.png", width = 4, height = 4, units = "in", res = 600)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0), cex.main = 0.9, cex.lab = 0.9, cex.axis = 0.9)
plot_support(fit_dpe)
dev.off()
