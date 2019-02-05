sep_sslg <- function(X, Y, 
                     lambdas = list(lambda1 = 1, lambda0 = seq(1, nrow(X), length = 10)), 
                     xis = list(xi1 = 0.01 * nrow(X), xi0 = seq(0.1 * nrow(X), nrow(X), length = 10)),
                     theta_hyper_params = c(1, ncol(X) * ncol(Y)), 
                     eta_hyper_params = c(1, ncol(Y)), 
                     diag_penalty = 1,
                     max_iter = 500,
                     eps = 1e-3,
                     obj_counter_max = 5, 
                     verbose = 0){
  n <- nrow(Y)
  q <- ncol(Y)
  p <- ncol(X)
  
  B_hat <- matrix(nrow = p, ncol = q)
  lambda1 <- lambdas[["lambda1"]]
  lambda0 <- lambdas[["lambda0"]]
  L <- length(lambda0)
  xi1 <- xis[["xi1"]]
  xi0 <- xis[["xi0"]]
  
  a_theta <- theta_hyper_params[1]
  b_theta <- theta_hyper_params[2]
  
  a_eta <- eta_hyper_params[1]
  b_eta <- eta_hyper_params[2]

  if(verbose == 1) print("Starting variable selection:")
  time_start <- Sys.time()
  for(k in 1:q){
    tmp <- SSLASSO(X, Y[,k], penalty = "adaptive", variance= "fixed", 
                   lambda1 = lambda1, lambda0 = lambda0, sigma = 1,
                   a = a_theta, b = b_theta, eps = eps, max.iter = max_iter)
    B_hat[,k] <- tmp$beta[,L]
  }
  if(verbose == 1) print(paste("  num B != 0 :", sum(B_hat != 0)))
  
  R <- Y - X %*% B_hat
  if(verbose == 1) print("Starting covariance selection:")
  #tmp <- gSSL(R, xis, eta_hyper_params, diag_penalty, control_params, verbose = verbose)
  tmp <- gSSL(R, xis, eta_hyper_params, diag_penalty, max_iter, eps, obj_counter_max,verbose)
  Omega_hat <- tmp$Omega
  time_end <- Sys.time()
  
  return(list(B = B_hat, Omega = Omega_hat, time = as.numeric(time_end - time_start)))
}