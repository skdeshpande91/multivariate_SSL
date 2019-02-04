\name{sep_sslg}
\alias{sep_sslg}
\title{Separate Variable and Covariance Selection in Multivariate Regression}
\description{
Consider the multivariate linear regression model Y = XB + E, where Y is an n x q matrix of outcomes, X is an n x p matrix of covariates, and the rows are i.i.d multivariate normal N(0, Omega^-1).
Rather than perform simultaneous variable and covariance selection, \code{sep_sslg} performs variable selection with SSL regularization for each outcome separately using the function \code{\link{SSLASSO}}.
Then it estimates a sparse precision matrix from the resulting matrix of residuals with SSL regularization using the function \code{\link{gSSL}}.
This procedure is described in more detail in Deshpande et al. (2018). 
Note: while this is much faster than \code{\link{mSSL_dpe}}, the final estimates of the matrix of covariate effects B and the residual precision matrix Omega are not MAP estimates. 
Indeed, the support recovery and estimation performance of \code{sep_sslg} is generally worse than that of \code{\link{mSSL_dpe}}. 
This function is provided for the sake of completeness and in the interest of reproducing examples in Deshpande et al. (2018). 
}
\usage{
sep_sslg(X, Y, lambdas = list(lambda1 = 1,lambda0 = seq(1, nrow(X), length = 10)), 
         xis = list(xi1 = 0.01 * nrow(X), xi0 = seq(0.1 * nrow(X),nrow(X), length = 10)),
         theta_hyper_params = c(1, ncol(X)*ncol(Y)), 
         eta_hyper_params = c(1, ncol(Y)), diag_penalty = 1,
         max_iter = 500,
         eps = 0.001, 
         verbose = 0)
}
\arguments{
  \item{X}{Matrix of covariates (n x p). Columns are centered and scaled to have norm sqrt(n) by default.}
  \item{Y}{Matrix of outcomes (n x q). Columns are centered by default.}
  \item{lambdas}{A list containing two elements:
    \itemize{
      \item lambda1 : Scale parameter of the slab distribution for the elements of B. Default value is 1. 
      \item lambda0 : Vector of length L, specifying the increasing sequence of spike distribution scale parameters for elements of B. Default is 10 evenly spaced values between 10 and n. 
    }
  }
  \item{xis}{A list containing two elements:
    \itemize{
      \item xi1 : Scale parameter of the slab distribution for the off-diagonal elements of Omega. Default value is 0.01 * n.
      \item xi0 : Vector of length L, specifying the increasing sequence of spike distribution scale parameters for the off-diagonal elements of Omega. Default is 10 evenly spaced values between 0.1 * n and n. 
    }
  }
  \item{theta_hyper_params}{Vector of length two specifying the hyper-parameters for the Beta prior on theta, the mixing weight of the spike-and-slab prior for elements of B. Default is a Beta(1, p*q) prior.}
  \item{eta_hyper_params}{Vector of length two specifying the hyper-parameters for the Beta prior on eta, the mixing weight of the spike-and-slab prior for off-diagonal elements of Omega. Default is a Beta(1, q) prior.}
  \item{diag_penalty}{Binary indicator whether or not the diagonal elements of Omega should be penalized (\code{diag_penalty = 1}). Default is 1.}
  \item{control_params}{A list containing the following elements used to monitor convergence:
    \itemize{
      \item max_iter : Maximum number of iterations for \code{\link{SSLASSO}} and \code{\link{gSSL}}. Default is 500.
      \item eps : Convergence tolerance for both \code{\link{SSLASSO}} and \code{\link{gSSL}}. Default is 0.001. 
      \item obj_counter_max : When the percentage increase in the objective function is less than \code{eps} for \code{obj_counter_max} consecutive iterations, the EM algorithm used in \code{\link{gSSL}} terminates. Default value is 5.
    }
  }
  \item{verbose}{Integer (0 or 1) that indicates whether or not to print the progress of the algorithm. Default is 0.}
}
\value{
  An list containing:
  \item{B}{The final estimate of covariate effects B.}
  \item{Omega}{The final estimate of the residual precision matrix Omega}
  \item{time}{Total time in seconds.}
 }
\references{Deshpande, S.K., Rockova, V., and George, E.I. (2018). "Simultaneous Variable and Covariance Selection with the Multivariate Spike-and-Slab LASSO". arXiv:1708.08911}
\author{Sameer Deshpande}
\seealso{\code{\link{mSSL_dcpe}}}
\examples{
library(SSLASSO)
library(MASS)
library(mSSL)
n <- 100
p <- 50
q <- 25
rho <- 0.9
 
set.seed(129)
SigmaX <- 0.7^(abs(outer(1:p, 1:p, "-")))
X <- mvrnorm(n, rep(0, times = p), SigmaX)
B.orig <- matrix(sample(c(runif(floor(p*q*0.2), -2,2), 
            rep(0, times = p*q - floor(p*q*0.2)))), 
            nrow = p, ncol = q)
XB <- X \%*\% B.orig
Sigma_orig <- rho^(abs(outer(1:q, 1:q, "-")))
Omega_orig <- solve(Sigma_orig)
Omega_orig[abs(Omega_orig) < 1e-3] <- 0
set.seed(130)
E <- mvrnorm(n, rep(0, times = q), Sigma_orig)
Y <- XB + E
fit_sep_sslg <- sep_sslg(X, Y)
}