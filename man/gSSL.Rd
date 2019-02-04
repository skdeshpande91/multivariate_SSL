\name{gSSL}
\alias{gSSL}
\title{Gaussian Graphical Modeling with the Spike-and-Slab LASSO}
\description{Estimates a sparse precision matrix from multivariate Gaussian data with spike-and-slab LASSO regularization. 
Specifically, spike-and-slab LASSO priors are placed on the off-diagonal elements of the precision matrix Omega.
We then use an EM algorithm, which reduces to solving a series of graphical LASSO problems with self-adaptive penalties.
We solve these problems using the QUIC algorithm of Hsieh et al. (2014). We have used code available in 
QUIC package available on CRAN. 
We note that this EM algorithm is a special case of the ECM algorithm implemented in \code{\link{mSSL_dpe}}
}
\usage{
gSSL(Y, 
    xis = list(xi1 = 0.01 * nrow(Y), xi0 = seq(0.1 * nrow(Y), nrow(Y), length = 10)), 
     eta_hyper_params = c(1, ncol(Y)), 
     diag_penalty = 1, 
     max_iter = 500,
     eps = 1e-3,
     obj_counter_max = 5,
     verbose = 0)
}
\arguments{
  \item{Y}{Data matrix (n x q). }
  \item{xis}{A list containing two elements:
    \itemize{
      \item xi1 : Scale parameter of the slab distribution for the off-diagonal elements of Omega. Default is 0.01 * n.
      \item xi0 : Vector of length L, specifying the increasing sequence of spike distribution scale parameters for the off-diagonal elements of Omega. Default is 10 equally spaced values between 0.1 * n and n.
    }
  }
  \item{eta_hyper_params}{Vector of length two specifying the hyper-parameters for the Beta prior on eta, the mixing weight of the spike-and-slab prior for off-diagonal elements of Omega. By default the prior is Beta(1,q).}
  \item{diag_penalty}{Binary indicator whether or not the diagonal elements of Omega should be penalized (\code{diag_penalty = 1}). Default is 1.}
  \item{max_iter}{Maximum number of iterations for the main ECM algorithm. Default value is 500.}
  \item{eps}{Convergence tolerance for the EM algorithm. If the percentage change in every element of B and Omega is less than \code{eps}, the algorithm terminates. Default value is 0.001.}
  \item{obj_counter_max}{When the percentage increase in the objective function is less than \code{eps} for \code{obj_counter_max} iterations, the EM algorithm terminates. Default value is 5.}
  \item{verbose}{Binary indicator whether or not to print the progress of the algorithm. Default is 0.}
}
\value{
  An list containing:
  \item{Omega}{The estimated residual precision matrix.}
  \item{Sigma}{The estimated residual covariance matrix.}
  \item{eta}{The estimate of the spike-and-slab mixing weight eta.}
  \item{Omega_path}{Array of MAP estimates of Omega corresponding to each value of xi0.}
  \item{Sigma_path}{Array of MAP estimates of Sigma corresponding to each value of xi0.}
  \item{eta_path}{Vector of MAP estimates of eta corresponding to each value of xi0.}
 }
\references{
  \itemize{
    \item{Deshpande, S.K., Rockova, V., and George, E.I. (2018). "Simultaneous Variable and Covariance Selection with the Multivariate Spike-and-Slab LASSO". arXiv:1708.08911}

    \item{Hsieh, C-J, Sustick, M.A., Dhillon, I.S., and Ravikumar, P. (2011). "Sparse Inverse Covariance Matrix Estimation Using Quadratic Approximation." Advances in Neural Information Processing Systems 24, 2330--2338.}
  }
}
\author{Sameer Deshpande}
\examples{
library(SSLASSO)
library(MASS)
n <- 100
p <- 50
q <- 25
rho <- 0.9
 
set.seed(130)
Sigma_orig <- rho^(abs(outer(1:q, 1:q, "-")))
Omega_orig <- solve(Sigma_orig)
Omega_orig[abs(Omega_orig) < 1e-3] <- 0
E <- mvrnorm(n, rep(0, times = q), Sigma_orig)

fit <- gSSL(E)
}