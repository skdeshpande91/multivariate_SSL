# Script used for simulations
# These were run on a high-performance computing environment as a task array, which passed in two arguments
# sim_number corresponds to the simulation numbering in the main paper. 
# The 100 simulation replicates were divided into 10 batches of 10 datasets.

library(MASS)
library(R.utils)
library(MRCE)
library(lpSolve)
library(SSLASSO)
library(mSSL)
source("scripts/load_capme.R")
source("scripts/new_mrce_cv.R")
source("scripts/sep_lasso_glasso.R")

args <- commandArgs(TRUE)
sim_number <- as.numeric(args[1])
batch <- as.numeric(args[2])

if(sim_number == 1){
  n <- 400
  p <- 500
  q <- 25
  q <- 0.9
} else if(sim_number == 2){
  n <- 400
  p <- 500
  q <- 25
  rho <- 0.7  
} else if(sim_number == 3){
  n <- 400
  p <- 500
  q <- 25
  rho <- 0.5  
} else if(sim_number == 4){
  n <- 400
  p <- 500
  q <- 25
  rho <- 0
} else if(sim_number == 5){
  n <- 100
  p <- 50
  q <- 25
  rho <- 0.9
} else if(sim_number == 6){
  n <- 100
  p <- 50
  q <- 25
  rho <- 0.7
} else if(sim_number == 7){
  n <- 100
  p <- 50
  q <- 25
  rho <- 0.5
} else if(sim_number == 8){
  n <- 100
  p <- 50
  q <- 25
  rho <- 0
}

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
# Define lists to store results for each method considered

mSSL_dpe_perf <- list()
mSSL_dpe_perf[["B"]] <-  matrix(nrow = reps, ncol = 12, dimnames = list(c(),c("TP", "TN", "FP", "FN", "SEN", "SPE", "PREC", "ACC",  "F1", "MCC", "MSE", "TIME")))
mSSL_dpe_perf[["Omega"]] <- matrix(nrow = reps, ncol = 12, dimnames = list(c(), c("TP", "TN", "FP", "FN", "SEN", "SPE", "PREC", "ACC",  "F1", "MCC", "FROB", "TIME")))

mSSL_dcpe_perf <- mSSL_dpe_perf
mrce_perf <- mSSL_dpe_perf
capme_perf <- mSSL_dpe_perf
sep_lg_perf <- mSSL_dpe_perf
sep_sslg_perf <- mSSL_dpe_perf

seed_seq <- c(129, 603, 211, 724, 927, 1123, 1991, 691, 7591, 3219)

# Set lambdas for MRCE
mrce_lam1_vec <- 10^seq(-3,0, length = 10)
mrce_lam2_vec <- 10^seq(-2,0,length = 10)

for(r in 1:reps){
  SEED <- seed_seq[batch] + r
  set.seed(SEED)
  E <- mvrnorm(n, rep(0, times = q), Sigma_orig)
  Y <- XB + E
  print(paste("starting r = ", r, "at", Sys.time()))
  time_mSSL_dpe <- tryCatch({withTimeout(system.time(fit_mSSL_dpe <- mSSL_dpe(X, Y)),timeout = 14400, onTimeout = "warning")}, error = function(e) NULL)
  print(paste("    Finished mEMVS-DPE at", Sys.time()))
  time_mSSL_dcpe <- tryCatch({withTimeout(system.time(fit_mSSL_dcpe <- mSSL_dcpe(X, Y)),timeout = 7200, onTimeout = "warning")}, error = function(e) NULL)
  print(paste("    Finishd meMVS-DCPE at", Sys.time()))

  time_mrce <- tryCatch({withTimeout(system.time(fit_mrce <- new.mrce.cv(X, Y, lam.vec.1 = mrce_lam1_vec, lam.vec.2 = mrce_lam2_vec)), timeout = 9000, onTimeout = "warning")}, error = function(e) NULL)
  print(paste("    Finished MRCE at", Sys.time()))
  time_capme <- tryCatch({withTimeout(system.time(fit_capme <- my.capme(X, Y, linsolver.Gamma = "simplex", linsolver.Omega = "simplex")),timeout = 9000, onTimeout = "warning")}, error = function(e) NULL)
  print(paste("    Finished CAPME at", Sys.time()))
  time_sep_lg <- tryCatch({withTimeout(system.time(fit_sep_lg <- sep.lasso.glasso(Y,X)), timeout = 3600, onTimeout = "warning")}, error = function(e) NULL)
  print(paste("    Finished Sep Lasso + Glasso at", Sys.time()))
  time_sep_sslg <- tryCatch({withTimeout(system.time(fit_sep_sslg <- sep_sslg(X, Y)), timeout = 3600, onTimeout = "warning")}, error = function(e) NULL)
  print(paste("    Finished Sep SSL + SSGlasso at", Sys.time()))
  
  if(!is.null(time_mSSL_dpe)){
    #print("non-null answer for DPE")
    B_hat_mSSL_dpe <- fit_mSSL_dpe$B
    Omega_hat_mSSL_dpe <- fit_mSSL_dpe$Omega
    mSSL_dpe_perf[["B"]][r,] <- c(error_B(B_hat_mSSL_dpe, B_orig), time_mSSL_dpe["elapsed"])
    mSSL_dpe_perf[["Omega"]][r,] <- c(error_Omega(Omega_hat_mSSL_dpe, Omega_orig), time_mSSL_dpe["elapsed"])
  }
  if(!is.null(time_mSSL_dcpe)){
    B_hat_mSSL_dcpe <- fit_mSSL_dcpe$B
    Omega_hat_mSSL.dcpe <- fit_mSSL_dcpe$Omega
    mSSL_dcpe_perf[["B"]][r,] <- c(error_B(B_hat_mSSL_dcpe, B_orig), time_mSSL_dcpe["elapsed"])
    mSSL_dcpe_perf[["Omega"]][r,] <- c(error_Omega(Omega_hat_mSSL_dcpe, Omega_orig),time_mSSL_dcpe["elapsed"])
  }
  if(!is.null(time_mrce)){
    B_hat_mrce <- fit_mrce$Bhat
    Omega_hat_mrce <- fit_mrce$omega
    mrce_perf[["B"]][r,] <- c(error_B(B_hat_mrce, B_orig), time_mrce["elapsed"])
    mrce_perf[["Omega"]][r,] <- c(error_Omega(Omega_hat_mrce, Omega_orig), time_mrce["elapsed"])
  }
  if(!is.null(time_capme)){
    B_hat_capme <- fit_capme$B
    Omega_hat_capme <- fit_capme$Omega
    capme_perf[["B"]][r,] <- c(error_B(B_hat_capme, B_orig), time_capme["elapsed"])
    capme_perf[["Omega"]][r,] <- c(error_Omega(Omega_hat_capme, Omega_orig), time_capme["elapsed"])
  }
  if(!is.null(time_sep_lg)){
    B_hat_sep_lg <- fit_sep_lg$B
    Omega_hat_sep_lg <- fit_sep_lg$Omega
    sep_lg_perf[["B"]][r,] <- c(error_B(B_hat_sep_lg, B_orig), time_sep_lg["elapsed"])
    sep_lg_perf[["Omega"]][r,] <- c(error_Omega(Omega_hat_sep_lg, Omega_orig), time_sep_lg["elapsed"])
  }
  if(!is.null(time_sep_sslg)){
    B_hat_sep_sslg <- fit_sep_sslg$B
    Omega_hat_sep_sslg <- fit_sep_sslg$Omega
    sep_sslg_perf[["B"]][r,] <- c(error_B(B_hat_sep_sslg, B_orig), time_sep_sslg["elapsed"])
    sep_sslg_perf[["Omega"]][r,] <- c(error_Omega(Omega_hat_sep_sslg, Omega_orig), time_sep_sslg["elapsed"])
  }
  
  print(paste("Finished r =", r, "at", Sys.time()))
}
assign(paste0("mSSL_dpe_sim", sim_number, ".", batch), mSSL_dpe_perf)
assign(paste0("mSSL_dcpe_sim", sim_number, ".", batch), mSSL_dcpe_perf)
assign(paste0("mrce_sim", sim_number, ".", batch), mrce_perf)
assign(paste0("capme_sim", sim_number, ".", batch), capme_perf)
assign(paste0("sep_lg_sim", sim_number, ".", batch), sep_lg_perf)
assign(paste0("sep_sslg_sim", sim_number, ".", batch), sep_sslg_perf)

save_list <- paste0(c("mSSL_dpe_sim", "mSSL_dcpe_sim", "mrce_sim", "capme_sim", "sep_lg_sim", "sep_sslg_sim"), sim_number, ".", batch)
save_list <- c(save_list, "B_orig", "Omega_orig", "n", "p", "q")

save(list = save_list, file = paste0("sim_results/sim", sim_number, "_", batch, ".RData"))

