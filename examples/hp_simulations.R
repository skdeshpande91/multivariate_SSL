# Simulations for assessing hyperparameter sensitivity
library(MASS)
library(SSLASSO)
library(mSSL)

args <- commandArgs(TRUE)
batch <- as.numeric(args[1])


n <- 100
p <- 50
q <- 25
rho <- 0.9

# Create X matrices
set.seed(129)
SigmaX <- 0.7^(abs(outer(1:p, 1:p, "-")))
X <- mvrnorm(n, rep(0, times = p), SigmaX)

B_orig <- matrix(sample(c(runif(floor(p*q*0.2), -2,2), rep(0, times = p*q - floor(p*q*0.2)))), nrow = p, ncol = q)
XB <- X %*% B_orig

Sigma_orig <- rho^(abs(outer(1:q, 1:q, "-")))
Omega_orig <- solve(Sigma_orig)
Omega_orig[abs(Omega_orig) < 1e-3] <- 0

reps <- 10
seed_seq <- c(129, 603, 211, 724, 927, 1123, 1991, 691, 7591, 3219)

#######################
# lambda1 = 1
# max(lambda0) = max(abs(t(X) %*% Y))
# spacing: on absolute scale or on log-scale
# min(lambda0) = 10 (just for separation. totally arbitrary but we don't really care about solution here anyway)
# 
# min(xi0) = max(abs(t(Y)%*%Y/n))/100 ... this is similar to what is done in the glasso package
# max(xi0) = max(abs(t(Y) %*% Y)/n)/10 
# xi1 = max(abs(t(Y) %*% Y)/n)/1000 ... again, just starting xi0 and xi1 far apart
#######################


# Consider several different schemes
# Scheme 1: Original (as in paper)
# Scheme 2: Original lambdas and xis, (b_theta = p, b_eta = q)
# Scheme 3: Original lambdas and xis, (b_theta = 1, b_eta = 1)
# Scheme 4: Original penalties, but on log-scale, original b parameters
# Scheme 5: Original penalties, on log-scale, (b_theta = p, b_eta = q)
# Scheme 6: Original penalties, on log-scale, (b_theta = 1, b_eta = 1)
# Scheme 7: new penalties, same b_eta and b_theta's as in paper
# Scheme 8: new penalties, (b_theta = p, b_eta = q)
# Scheme 9: new penalties, (b_theta = 1, b_eta = 1)
# Scheme 10: new penalties, on log-scale, (b_theta = pq, b_eta = q)
# Scheme 11: new penalties, on log-scale, (b_theta = p, b_eta = q)
# Scheme 12: new penalties, on log-scale, (b_theta = 1, b_eta = 1)

fit1_dpe_perf <- list()
fit1_dpe_perf[["B"]] <-  matrix(nrow = reps, ncol = 12, dimnames = list(c(),c("TP", "TN", "FP", "FN", "SEN", "SPE", "PREC", "ACC",  "F1", "MCC", "MSE", "TIME")))
fit1_dpe_perf[["Omega"]] <- matrix(nrow = reps, ncol = 12, dimnames = list(c(), c("TP", "TN", "FP", "FN", "SEN", "SPE", "PREC", "ACC",  "F1", "MCC", "FROB", "TIME")))

fit1_dcpe_perf <- fit1_dpe_perf

fit2_dpe_perf <- fit1_dpe_perf
fit2_dcpe_perf <- fit1_dpe_perf
fit3_dpe_perf <- fit1_dpe_perf
fit3_dcpe_perf <- fit1_dpe_perf
fit4_dpe_perf <- fit1_dpe_perf
fit4_dcpe_perf <- fit1_dpe_perf
fit5_dpe_perf <- fit1_dpe_perf
fit5_dcpe_perf <- fit1_dpe_perf
fit6_dpe_perf <- fit1_dpe_perf
fit6_dcpe_perf <- fit1_dpe_perf
fit7_dpe_perf <- fit1_dpe_perf
fit7_dcpe_perf <- fit1_dpe_perf
fit8_dpe_perf <- fit1_dpe_perf
fit8_dcpe_perf <- fit1_dpe_perf
fit9_dpe_perf <- fit1_dpe_perf
fit9_dcpe_perf <- fit1_dpe_perf
fit10_dpe_perf <- fit1_dpe_perf
fit10_dcpe_perf <- fit1_dpe_perf
fit11_dpe_perf <- fit1_dpe_perf
fit11_dcpe_perf <- fit1_dpe_perf
fit12_dpe_perf <- fit1_dpe_perf
fit12_dcpe_perf <- fit1_dpe_perf

max_lam0_list <- rep(NA, times = reps)
max_xi0_list <- rep(NA, times = reps)
L <- 10

for(r in 1:reps){
  SEED <- seed_seq[batch] + r
  set.seed(SEED)
  E <- mvrnorm(n, rep(0, times = q), Sigma_orig)
  Y <- XB + E
  
  Y_cent - scale(Y, center = TRUE, scale = TRUE)
  max_lam0 <- max(abs(t(X) %*% Y))
  max_xi0 <- max(abs(n^-1 * t(Y_cent) %*% Y_cent))
  
  max_lam0_list[r] <- max_lam0
  max_xi0_list[r] <- max_xi0
  
  print(paste("starting r = ", r, "at", Sys.time()))
  # Scheme 1
  lambdas <- list(lambda1 = 1, lambda0 = seq(10, n ,length = L))
  xis <- list(xi1 = 0.01 * n, xi0 = n*seq(0.1, 1, length = L))
  theta_hyper_parameters <- c(1,p*q)
  eta_hyper_parameters <-  c(1, q)
  time_fit1_dpe <- tryCatch({withTimeout(system.time(fit1_dpe <- mSSL_dpe(X, Y, lambdas, xis, theta_hyper_params, eta_hyper_params)), timeout = 3600, onTimeout = "warning")}, error = function(e) NULL, warning = function(e) NULL)
  print(paste("    Finished Scheme 1 DPE at", Sys.time()))
  time_fit1_dcpe <- tryCatch({withTimeout(system.time(fit1_dcpe <- mSSL_dcpe(X, Y, lambdas, xis, theta_hyper_params, eta_hyper_params,1,0)), timeout = 3600, onTimeout = "warning")}, error = function(e) NULL, warning = function(e) NULL)  
  print(paste("    Finished Scheme1 DCPE at", Sys.time()))
  
  
  # Scheme 2
  lambdas <- list(lambda1 = 1, lambda0 = seq(10, n ,length = L))
  xis <- list(xi1 = 0.01 * n, xi0 = n*seq(0.1, 1, length = L))
  theta_hyper_params <- c(1, p)
  eta_hyper_params <- c(1, q)
  time_fit2_dpe <- tryCatch({withTimeout(system.time(fit2_dpe <- mSSL_dpe(X, Y, lambdas, xis, theta_hyper_params, eta_hyper_params)), timeout = 3600, onTimeout = "warning")}, error = function(e) NULL, warning = function(e) NULL)
  print(paste("    Finished Scheme 2 DPE at", Sys.time()))
  time_fit2_dcpe <- tryCatch({withTimeout(system.time(fit2_dcpe <- mSSL_dcpe(X, Y, lambdas, xis, theta_hyper_params, eta_hyper_params)), timeout = 3600, onTimeout = "warning")}, error = function(e) NULL, warning = function(e) NULL)  
  print(paste("    Finished Scheme 2 DCPE at", Sys.time()))
  
  # Scheme 3
  lambdas <- list(lambda1 = 1, lambda0 = seq(10, n ,length = L))
  xis <- list(xi1 = 0.01 * n, xi0 = n*seq(0.1, 1, length = L))
  theta_hyper_params <- c(1, 1)
  eta_hyper_params <- c(1, 1)
  time_fit3_dpe <- tryCatch({withTimeout(system.time(fit3_dpe <- mSSL_dpe(X, Y, lambdas, xis, theta_hyper_params, eta_hyper_params)), timeout = 3600, onTimeout = "warning")}, error = function(e) NULL, warning = function(e) NULL)
  print(paste("    Finished Scheme 3 DPE at", Sys.time()))
  time_fit3_dcpe <- tryCatch({withTimeout(system.time(fit3_dcpe <- mSSL_dcpe(X, Y, lambdas, xis, theta_hyper_params, eta_hyper_params)),timeout = 3600, onTimeout = "warning")}, error = function(e) NULL, warning = function(e) NULL)  
  print(paste("    Finished Scheme 3 DCPE at", Sys.time()))
  
  # Scheme 4
  lambdas <- list(lambda1 = 1, lambda0 = 10^seq(log(10, base = 10), log(n, base = 10) ,length = L))
  xis <- list(xi1 = 0.01 * n, xi0 = n*10^seq(log(0.1, base = 10), log(1, base = 10), length = L))
  theta_hyper_params <- c(1, p*q)
  eta_hyper_params <- c(1, q)
  time_fit4_dpe <- tryCatch({withTimeout(system.time(fit4_dpe <- mSSL_dpe(X, Y, lambdas, xis, theta_hyper_params, eta_hyper_params)), timeout = 3600, onTimeout = "warning")}, error = function(e) NULL, warning = function(e) NULL)
  print(paste("    Finished Scheme 4 DPE at", Sys.time()))
  time_fit4_dcpe <- tryCatch({withTimeout(system.time(fit4_dcpe <- mSSL_dcpe(X, Y, lambdas, xis, theta_hyper_params, eta_hyper_params)), timeout = 3600, onTimeout = "warning")}, error = function(e) NULL, warning = function(e) NULL)  
  print(paste("    Finished Scheme 4 DCPE at", Sys.time()))
  
  # Scheme 5
  lambdas <- list(lambda1 = 1, lambda0 = 10^seq(log(10, base = 10), log(n, base = 10) ,length = L))
  xis <- list(xi1 = 0.01 * n, xi0 = n*10^seq(log(0.1, base = 10), log(1, base = 10), length = L))
  theta_hyper_params <- c(1, p)
  eta_hyper_params <- c(1, q)
  time_fit5_dpe <- tryCatch({withTimeout(system.time(fit5_dpe <- mSSL_dpe(X, Y, lambdas, xis, theta_hyper_params, eta_hyper_params)), timeout = 3600, onTimeout = "warning")}, error = function(e) NULL, warning = function(e) NULL)
  print(paste("    Finished Scheme 5 DPE at", Sys.time()))
  time_fit5_dcpe <- tryCatch({withTimeout(system.time(fit5_dcpe <- mSSL_dcpe(X, Y, lambdas, xis, theta_hyper_params, eta_hyper_params)), timeout = 3600, onTimeout = "warning")}, error = function(e) NULL, warning = function(e) NULL)
  print(paste("    Finished Scheme 5 DCPE at", Sys.time()))  
  
  # Scheme 6
  lambdas <- list(lambda1 = 1, lambda0 = 10^seq(log(10, base = 10), log(n, base = 10) ,length = L))
  xis <- list(xi1 = 0.01 * n, xi0 = n*10^seq(log(0.1, base = 10), log(1, base = 10), length = L))
  theta_hyper_params <- c(1, 1)
  eta_hyper_params <- c(1, 1)
  time_fit6_dpe <- tryCatch({withTimeout(system.time(fit6_dpe <- mSSL_dpe(X, Y, lambdas, xis, theta_hyper_params, eta_hyper_params)),timeout = 3600, onTimeout = "warning")}, error = function(e) NULL, warning = function(e) NULL)
  print(paste("    Finished Scheme 6 DPE at", Sys.time()))
  time_fit6_dcpe <- tryCatch({withTimeout(system.time(fit6_dcpe <- mSSL_dcpe(X, Y, lambdas, xis, theta_hyper_params, eta_hyper_params)), timeout = 3600, onTimeout = "warning")}, error = function(e) NULL, warning = function(e) NULL)
  print(paste("    Finished Scheme 6 DCPE at", Sys.time()))
  
  
  # Scheme 7
  
  lambdas <- list(lambda1 = 1, lambda0 = seq(10, max_lam0 ,length = L))
  xis <- list(xi1 = n*max_xi0/1000, xi0 = n*seq(max_xi0/100, max_xi0/10, length = L))
  theta_hyper_params <- c(1, p*q)
  eta_hyper_params <- c(1, q)
  time_fit7_dpe <- tryCatch({withTimeout(system.time(fit7_dpe <- mSSL_dpe(X, Y, lambdas, xis, theta_hyper_params, eta_hyper_params)), timeout = 3600, onTimeout = "warning")}, error = function(e) NULL, warning = function(e) NULL)
  print(paste("    Finished Scheme 7 DPE at", Sys.time()))
  time_fit7_dcpe <- tryCatch({withTimeout(system.time(fit7_dcpe <- mSSL_dcpe(X, Y, lambdas, xis, theta_hyper_params, eta_hyper_params,1,0)), timeout = 3600, onTimeout = "warning")}, error = function(e) NULL, warning = function(e) NULL) 
  print(paste("    Finished Scheme 7 DCPE at", Sys.time()))
  
  # Scheme 8
  lambdas <- list(lambda1 = 1, lambda0 = seq(10, max_lam0 ,length = L))
  xis <- list(xi1 = n*max_xi0/1000, xi0 = n*seq(max_xi0/100, max_xi0/10, length = L))
  theta_hyper_params <- c(1, p)
  eta_hyper_params <- c(1, q)
  time_fit8_dpe <- tryCatch({withTimeout(system.time(fit8_dpe <- mSSL_dpe(X, Y, lambdas, xis, theta_hyper_params, eta_hyper_params)), timeout = 3600, onTimeout = "warning")}, error = function(e) NULL, warning = function(e) NULL)
  print(paste("    Finished Scheme 8 DPE at", Sys.time()))
  time_fit8_dcpe <- tryCatch({withTimeout(system.time(fit8_dcpe <- mSSL_dcpe(X, Y, lambdas, xis, theta_hyper_params, eta_hyper_params)), timeout = 3600, onTimeout = "warning")}, error = function(e) NULL, warning = function(e) NULL) 
  print(paste("    Finished Scheme 8 DCPE at", Sys.time()))
  
  
  # Scheme 9
  lambdas <- list(lambda1 = 1, lambda0 = seq(10, max_lam0 ,length = L))
  xis <- list(xi1 = n*max_xi0/1000, xi0 = n*seq(max_xi0/100, max_xi0/10, length = L))
  theta_hyper_params <- c(1, 1)
  eta_hyper_params <- c(1, 1)
  time_fit9_dpe <- tryCatch({withTimeout(system.time(fit9_dpe <- mSSL_dpe(X, Y, lambdas, xis, theta_hyper_params, eta_hyper_params,1,0)), timeout = 3600, onTimeout = "warning")}, error = function(e) NULL, warning = function(e) NULL)
  print(paste("    Finished Scheme 9 DPE at", Sys.time()))
  time_fit9_dcpe <- tryCatch({withTimeout(system.time(fit9_dcpe <- mSSL_dcpe(X, Y, lambdas, xis, theta_hyper_params, eta_hyper_params,1,0)), timeout = 3600, onTimeout = "warning")}, error = function(e) NULL, warning = function(e) NULL)
  print(paste("    Finished Scheme 9 DCPE at", Sys.time()))
  
  # Scheme 10
  lambdas <- list(lambda1 = 1, lambda0 = 10^seq(log(10, base = 10), log(max_lam0, base = 10),length = L))
  xis <- list(xi1 = n*max_xi0/1000, xi0 = n*10^seq(log(max_xi0/100, base = 10), log(max_xi0/10, base = 10), length = L))
  theta_hyper_params <- c(1, p*q)
  eta_hyper_params <- c(1, q)
  time_fit10_dpe <- tryCatch({withTimeout(system.time(fit10_dpe <- mSSL_dpe(X, Y, lambdas, xis, theta_hyper_params, eta_hyper_params)), timeout = 3600, onTimeout = "warning")}, error = function(e) NULL, warning = function(e) NULL)
  print(paste("    Finished Scheme 10 DPE at", Sys.time()))
  time_fit10_dcpe <- tryCatch({withTimeout(system.time(fit10_dcpe <- mSSL_dcpe(X, Y, lambdas, xis, theta_hyper_params, eta_hyper_params)), timeout = 3600, onTimeout = "warning")}, error = function(e) NULL, warning = function(e) NULL) 
  print(paste("    Finished Scheme 10 DCPE at", Sys.time()))
  
  # Scheme 11
  lambdas <- list(lambda1 = 1, lambda0 = 10^seq(log(10, base = 10), log(max_lam0, base = 10),length = L))
  xis <- list(xi1 = n*max_xi0/1000, xi0 = n*10^seq(log(max_xi0/100, base = 10), log(max_xi0/10, base = 10), length = L))
  theta_hyper_params <- c(1, p)
  eta_hyper_params <- c(1, q)
  time_fit11_dpe <- tryCatch({withTimeout(system.time(fit11_dpe <- mSSL_dpe(X, Y, lambdas, xis, theta_hyper_params, eta_hyper_params)), timeout = 3600, onTimeout = "warning")}, error = function(e) NULL, warning = function(e) NULL)
  print(paste("    Finished Scheme 11 DPE at", Sys.time()))
  time_fit11_dcpe <- tryCatch({withTimeout(system.time(fit11_dcpe <- mSSL_dcpe(X, Y, lambdas, xis, theta_hyper_params, eta_hyper_params)), timeout = 3600, onTimeout = "warning")}, error = function(e) NULL, warning = function(e) NULL) 
  print(paste("    Finished Scheme 11 DCPE at", Sys.time()))  
  
  # Scheme 12
  lambdas <- list(lambda1 = 1, lambda0 = 10^seq(log(10, base = 10), log(max_lam0, base = 10),length = L))
  xis <- list(xi1 = n*max_xi0/1000, xi0 = n*10^seq(log(max_xi0/100, base = 10), log(max_xi0/10, base = 10), length = L))
  theta_hyper_params <- c(1, 1)
  eta_hyper_params <- c(1, 1)
  time_fit12_dpe <- tryCatch({withTimeout(system.time(fit12_dpe <- mSSL_dpe(X, Y, lambdas, xis, theta_hyper_params, eta_hyper_params)), timeout = 3600, onTimeout = "warning")}, error = function(e) NULL, warning = function(e) NULL)
  print(paste("    Finished Scheme 12 DPE at", Sys.time()))
  time_fit12_dcpe <- tryCatch({withTimeout(system.time(fit12_dcpe <- mSSL_dcpe(X, Y, lambdas, xis, theta_hyper_params, eta_hyper_params)),  timeout = 3600, onTimeout = "warning")}, error = function(e) NULL, warning = function(e) NULL)
  print(paste("    Finished Scheme 12 DCPE at", Sys.time()))
  
  if(!is.null(time_fit1_dpe)){
    fit1_dpe_perf[["B"]][r,] <- c(error_B(fit1_dpe$B, B_orig), time_fit1_dpe["elapsed"])
    fit1_dpe_perf[["Omega"]][r,] <- c(error_Omega(fit1_dpe$Omega, Omega_orig), time_fit1_dpe["elapsed"])
  }
  if(!is.null(time_fit1_dcpe)){
    fit1_dcpe_perf[["B"]][r,] <- c(error_B(fit1_dcpe$B, B_orig), time_fit1_dcpe["elapsed"])
    fit1_dcpe_perf[["Omega"]][r,] <- c(error_Omega(fit1_dcpe$Omega, Omega_orig), time_fit1_dcpe["elapsed"])
  }
  if(!is.null(time_fit2_dpe)){
    fit2_dpe_perf[["B"]][r,] <- c(error_B(fit2_dpe$B, B_orig), time_fit2_dpe["elapsed"])
    fit2_dpe_perf[["Omega"]][r,] <- c(error_Omega(fit2_dpe$Omega, Omega_orig), time_fit2_dpe["elapsed"])
  }
  if(!is.null(time_fit2_dcpe)){
    fit2_dcpe_perf[["B"]][r,] <- c(error_B(fit2_dcpe$B, B_orig), time_fit2_dcpe["elapsed"])
    fit2_dcpe_perf[["Omega"]][r,] <- c(error_Omega(fit2_dcpe$Omega, Omega_orig), time_fit2_dcpe["elapsed"])
  }
  
  if(!is.null(time_fit3_dpe)){
    fit3_dpe_perf[["B"]][r,] <- c(error_B(fit3_dpe$B, B_orig), time_fit3_dpe["elapsed"])
    fit3_dpe_perf[["Omega"]][r,] <- c(error_Omega(fit3_dpe$Omega, Omega_orig), time_fit3_dpe["elapsed"])
  }
  if(!is.null(time_fit3_dcpe)){
    fit3_dcpe_perf[["B"]][r,] <- c(error_B(fit3_dcpe$B, B_orig), time_fit3_dcpe["elapsed"])
    fit3_dcpe_perf[["Omega"]][r,] <- c(error_Omega(fit3_dcpe$Omega, Omega_orig), time_fit3_dcpe["elapsed"])
  }
  
  if(!is.null(time_fit4_dpe)){
    fit4_dpe_perf[["B"]][r,] <- c(error_B(fit4_dpe$B, B_orig), time_fit4_dpe["elapsed"])
    fit4_dpe_perf[["Omega"]][r,] <- c(error_Omega(fit4_dpe$Omega, Omega_orig), time_fit4_dpe["elapsed"])
  }
  if(!is.null(time_fit4_dcpe)){
    fit4_dcpe_perf[["B"]][r,] <- c(error_B(fit4_dcpe$B, B_orig), time_fit4_dcpe["elapsed"])
    fit4_dcpe_perf[["Omega"]][r,] <- c(error_Omega(fit4_dcpe$Omega, Omega_orig), time_fit4_dcpe["elapsed"])
  }
  
  if(!is.null(time_fit5_dpe)){
    fit5_dpe_perf[["B"]][r,] <- c(error_B(fit5_dpe$B, B_orig), time_fit5_dpe["elapsed"])
    fit5_dpe_perf[["Omega"]][r,] <- c(error_Omega(fit5_dpe$Omega, Omega_orig), time_fit5_dpe["elapsed"])
  }
  if(!is.null(time_fit5_dcpe)){
    fit5_dcpe_perf[["B"]][r,] <- c(error_B(fit5_dcpe$B, B_orig), time_fit5_dcpe["elapsed"])
    fit5_dcpe_perf[["Omega"]][r,] <- c(error_Omega(fit5_dcpe$Omega, Omega_orig), time_fit5_dcpe["elapsed"])
  }
  
  if(!is.null(time_fit6_dpe)){
    fit6_dpe_perf[["B"]][r,] <- c(error_B(fit6_dpe$B, B_orig), time_fit6_dpe["elapsed"])
    fit6_dpe_perf[["Omega"]][r,] <- c(error_Omega(fit6_dpe$Omega, Omega_orig), time_fit6_dpe["elapsed"])
  }
  if(!is.null(time_fit6_dcpe)){
    fit6_dcpe_perf[["B"]][r,] <- c(error_B(fit6_dcpe$B, B_orig), time_fit6_dcpe["elapsed"])
    fit6_dcpe_perf[["Omega"]][r,] <- c(error_Omega(fit6_dcpe$Omega, Omega_orig), time_fit6_dcpe["elapsed"])
  }
  
  if(!is.null(time_fit7_dpe)){
    fit7_dpe_perf[["B"]][r,] <- c(error_B(fit7_dpe$B, B_orig), time_fit7_dpe["elapsed"])
    fit7_dpe_perf[["Omega"]][r,] <- c(error_Omega(fit7_dpe$Omega, Omega_orig), time_fit7_dpe["elapsed"])
  }
  if(!is.null(time_fit7_dcpe)){
    fit7_dcpe_perf[["B"]][r,] <- c(error_B(fit7_dcpe$B, B_orig), time_fit7_dcpe["elapsed"])
    fit7_dcpe_perf[["Omega"]][r,] <- c(error_Omega(fit7_dcpe$Omega, Omega_orig), time_fit7_dcpe["elapsed"])
  }
  
  if(!is.null(time_fit8_dpe)){
    fit8_dpe_perf[["B"]][r,] <- c(error_B(fit8_dpe$B, B_orig), time_fit8_dpe["elapsed"])
    fit8_dpe_perf[["Omega"]][r,] <- c(error_Omega(fit8_dpe$Omega, Omega_orig), time_fit8_dpe["elapsed"])
  }
  if(!is.null(time_fit8_dcpe)){
    fit8_dcpe_perf[["B"]][r,] <- c(error_B(fit8_dcpe$B, B_orig), time_fit8_dcpe["elapsed"])
    fit8_dcpe_perf[["Omega"]][r,] <- c(error_Omega(fit8_dcpe$Omega, Omega_orig), time_fit8_dcpe["elapsed"])
  }
  
  if(!is.null(time_fit9_dpe)){
    fit9_dpe_perf[["B"]][r,] <- c(error_B(fit9_dpe$B, B_orig), time_fit9_dpe["elapsed"])
    fit9_dpe_perf[["Omega"]][r,] <- c(error_Omega(fit9_dpe$Omega, Omega_orig), time_fit9_dpe["elapsed"])
  }  
  
  if(!is.null(time_fit9_dcpe)){
    fit9_dcpe_perf[["B"]][r,] <- c(error_B(fit9_dcpe$B, B_orig), time_fit9_dcpe["elapsed"])
    fit9_dcpe_perf[["Omega"]][r,] <- c(error_Omega(fit9_dcpe$Omega, Omega_orig), time_fit9_dcpe["elapsed"])
  }
  if(!is.null(time_fit10_dpe)){
    fit10_dpe_perf[["B"]][r,] <- c(error_B(fit10_dpe$B, B_orig), time_fit10_dpe["elapsed"])
    fit10_dpe_perf[["Omega"]][r,] <- c(error_Omega(fit10_dpe$Omega, Omega_orig), time_fit10_dpe["elapsed"])
  }
  if(!is.null(time_fit10_dcpe)){
    fit10_dcpe_perf[["B"]][r,] <- c(error_B(fit10_dcpe$B, B_orig), time_fit10_dcpe["elapsed"])
    fit10_dcpe_perf[["Omega"]][r,] <- c(error_Omega(fit10_dcpe$Omega, Omega_orig), time_fit10_dcpe["elapsed"])
  }
  if(!is.null(time_fit11_dpe)){
    fit11_dpe_perf[["B"]][r,] <- c(error_B(fit11_dpe$B, B_orig), time_fit11_dpe["elapsed"])
    fit11_dpe_perf[["Omega"]][r,] <- c(error_Omega(fit11_dpe$Omega, Omega_orig), time_fit11_dpe["elapsed"])
  }
  if(!is.null(time_fit11_dcpe)){
    fit11_dcpe_perf[["B"]][r,] <- c(error_B(fit11_dcpe$B, B_orig), time_fit11_dcpe["elapsed"])
    fit11_dcpe_perf[["Omega"]][r,] <- c(error_Omega(fit11_dcpe$Omega, Omega_orig), time_fit11_dcpe["elapsed"])
  }
  if(!is.null(time_fit12_dpe)){
    fit12_dpe_perf[["B"]][r,] <- c(error_B(fit12_dpe$B, B_orig), time_fit12_dpe["elapsed"])
    fit12_dpe_perf[["Omega"]][r,] <- c(error_Omega(fit12_dpe$Omega, Omega_orig), time_fit12_dpe["elapsed"])
  }
  if(!is.null(time_fit12_dcpe)){
    fit12_dcpe_perf[["B"]][r,] <- c(error_B(fit12_dcpe$B, B_orig), time_fit12_dcpe["elapsed"])
    fit12_dcpe_perf[["Omega"]][r,] <- c(error_Omega(fit12_dcpe$Omega, Omega_orig), time_fit12_dcpe["elapsed"])
  }
  print(paste("Finished r =", r, "at", Sys.time()))
}

assign(paste0("fit1_dpe_", batch), fit1_dpe_perf)
assign(paste0("fit1_dcpe_", batch), fit1_dcpe_perf)
assign(paste0("fit2_dpe_", batch), fit2_dpe_perf)
assign(paste0("fit2_dcpe_", batch), fit2_dcpe_perf)
assign(paste0("fit3_dpe_", batch), fit3_dpe_perf)
assign(paste0("fit3_dcpe_", batch), fit3_dcpe_perf)
assign(paste0("fit4_dpe_", batch), fit4_dpe_perf)
assign(paste0("fit4_dcpe_", batch), fit4_dcpe_perf)
assign(paste0("fit5_dpe_", batch), fit5_dpe_perf)
assign(paste0("fit5_dcpe_", batch), fit5_dcpe_perf)
assign(paste0("fit6_dpe_", batch), fit6_dpe_perf)
assign(paste0("fit6_dcpe_", batch), fit6_dcpe_perf)
assign(paste0("fit7_dpe_", batch), fit7_dpe_perf)
assign(paste0("fit7_dcpe_", batch), fit7_dcpe_perf)
assign(paste0("fit8_dpe_", batch), fit8_dpe_perf)
assign(paste0("fit8_dcpe_", batch), fit8_dcpe_perf)
assign(paste0("fit9_dpe_", batch), fit9_dpe_perf)
assign(paste0("fit9_dcpe_", batch), fit9_dcpe_perf)
assign(paste0("fit10_dpe_", batch), fit10_dpe_perf)
assign(paste0("fit10_dcpe_", batch), fit10_dcpe_perf)
assign(paste0("fit11_dpe_", batch), fit11_dpe_perf)
assign(paste0("fit11_dcpe_", batch), fit11_dcpe_perf)
assign(paste0("fit12_dpe_", batch), fit12_dpe_perf)
assign(paste0("fit12_dcpe_", batch), fit12_dcpe_perf)
assign(paste0("max_lam0_", batch), max_lam0_list)
assign(paste0("max_xi0_", batch), max_xi0_list)


save_list_dpe <- paste0(c("fit1_dpe_", "fit2_dpe_", "fit3_dpe_", "fit4_dpe_", "fit5_dpe_", "fit6_dpe_", "fit7_dpe_", "fit8_dpe_", "fit9_dpe_", "fit10_dpe_", "fit11_dpe_", "fit12_dpe_"), batch)
save_list_dcpe <- paste0(c("fit1_dcpe_", "fit2_dcpe_", "fit3_dcpe_", "fit4_dcpe_", "fit5_dcpe_", "fit6_dcpe_", "fit7_dcpe_", "fit8_dcpe_", "fit9_dcpe_", "fit10_dcpe_", "fit11_dcpe_", "fit12_dcpe_"), batch)
save_list <- c(save_list_dpe, save_list_dcpe, paste0(c("max_lam0_", "max_xi0_"), batch), "B_orig", "Omega_orig", "n", "p", "q")

save(list = save_list, file = paste0("sim_results/hp_", batch, ".RData"))

