\name{mSSL_dcpe}
\alias{mSSL_dcpe}
\title{Dynamic Conditional Posterior Exploration with the multivariate spike-and-slab LASSO}
\description{
Consider the multivariate linear regression model Y = XB + E, where Y is an n x q matrix of outcomes, X is an n x p matrix of covariates, and the rows are i.i.d multivariate normal N(0, Omega^-1).
Like in \code{\link{mSSL_dpe}}, we target sparse estimates of B and Omega using spike-and-slab LASSO regularization.
Since identifying the joint MAP estimate of (B,Omega) can be computationally intensive, \code{mSSL_dcpe} first identifies a sequence of conditional maximizers of B, keeping Omega
fixed at the identity, indexed by an increasing of sequence of spike penalties.
Then, holding the value of B fixed at the terminal estimate of this sequence, \code{mSSL_dcpe} identifies a sequence of conditional maximizers of Omega indexed by another increasing sequence of spike penalties.
Finally, it refines the terminal conditonal MAP estimates of B and Omega to the joint MAP corresponding to the largest spike penalties for both B and Omega.
See Deshpande et al. (2018) for more details.
Note: \code{mSSL_dcpe} is substantially faster than \code{\link{mSSL_dpe}}. However, it tends to have worse support recovery and estimation performance. 
It is recommended to run both \code{mSSL_dcpe} and \code{\link{mSSL_dpe}} and to compare the results.
}
\usage{
mSSL_dcpe(X, Y, lambdas = list(lambda1 = 1, lambda0 = seq(10, nrow(X), length = 10)),
         xis = list(xi1 = 0.01 * nrow(X), xi0 = seq(0.1 * nrow(X), nrow(X), length = 10)),
         theta_hyper_params = c(1, ncol(X) * ncol(Y)), 
         eta_hyper_params = c(1, ncol(Y)),
         diag_penalty = 1, 
         max_iter = 500,
         eps = 1e-3,
         verbose = 0)
}
\arguments{
  \item{X}{Matrix of covariates (n x p). Columns are centered and scaled to have norm sqrt(n) by default.}
  \item{Y}{Matrix of outcomes (n x q). Columns are centered by default.}
  \item{lambdas}{A list containing two elements: 
    \itemize{
      \item lambda1 : Scale parameter of the slab distribution for the elements of B. Default value is 1.
      \item lambda0 : Vector of length L, specifying the increasing sequence of spike distribution scale parameters for elements of B. Default is sequence of 10 evenly spaced values between 10 and n.
    }
  }
  \item{xis}{A list containing two elements:
    \itemize{
      \item xi1 : Scale parameter of the slab distribution for the off-diagonal elements of Omega. Default value is 0.01 * n. 
      \item xi0 : Vector of length L, specifying the increasing sequence of spike distribution scale parameters for the off-diagonal elements of Omega. Default is sequence of 10 evenly spaced values between 0.1 * n and n.
    }
  }
  \item{theta_hyper_params}{Vector of length two specifying the hyper-parameters for the Beta prior on theta, the mixing weight of the spike-and-slab prior for elements of B. Default is a Beta(1, p*q) prior.}
  \item{eta_hyper_params}{Vector of length two specifying the hyper-parameters for the Beta prior on eta, the mixing weight of the spike-and-slab prior for off-diagonal elements of Omega. Default is a Beta(1, q) prior.}
  \item{diag_penalty}{Integer (0 or 1) that indicates whether or not the diagonal elements of Omega should be penalty (\code{diag_penalty = 1}). Default is 1.}
  \item{max_iter}{Maximum number of iterations for the main ECM algorithm. Default value is 500.}
  \item{eps}{Convergence tolerance for the ECM algorithm. If the percentage change in every element of B and Omega is less than \code{eps}, the algorithm terminates. Default value is 0.001.}
  \item{verbose}{Integer (0 or 1) that indicates whether or not to print the progress of the algorithm. Default is 0.}
}
\value{
  An list containing:
  \item{alpha}{The final estimate of the intercept}
  \item{B}{The final estimate of covariate effects B.}
  \item{Omega}{The final estimate of the residual precision matrix Omega}
  \item{theta}{The final estimate of theta}
  \item{eta}{The final estimate of eta}
  \item{alpha_path}{}
  \item{B0_path}{A p x q x L array containing conditional MAP estimate of B, conditional on Omega being the identity matrix, for a sequence of posteriors}
  \item{Omega0_path}{A q x q x L array containing conditional MAP estimates of Omega, conditional on B}
  \item{time}{Total time in seconds.}
 }
\references{
  \itemize{
    \item{Deshpande, S.K., Rockova, V., and George, E.I. (2018). "Simultaneous Variable and Covariance Selection with the Multivariate Spike-and-Slab LASSO". arXiv:1708.08911}

    \item{Hsieh, C-J, Sustick, M.A., Dhillon, I.S., and Ravikumar, P. (2014). "QUIC: Quadratic Approximation for Sparse Inverse Covariance Estimation." Journal of Machine Learning Research. 15, 2911 -- 2947}
  }
}
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
B_orig <- matrix(sample(c(runif(floor(p*q*0.2), -2,2), 
          rep(0, times = p*q - floor(p*q*0.2)))), 
          nrow = p, ncol = q)
XB <- X \%*\% B_orig
Sigma_orig <- rho^(abs(outer(1:q, 1:q, "-")))
Omega_orig <- solve(Sigma_orig)
Omega_orig[abs(Omega_orig) < 1e-3] <- 0

set.seed(130)
E <- mvrnorm(n, rep(0, times = q), Sigma_orig)
Y <- XB + E
fit_dcpe <- mSSL_dcpe(X, Y)

}
