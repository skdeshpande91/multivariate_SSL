\name{mSSL_dpe}
\alias{mSSL_dpe}
\title{Dynamic Posterior Exploration with the multivariate spike-and-slab LASSO}
\description{Consider the multivariate linear regression model Y = XB + E, where Y is an n x q matrix of outcomes, X is an n x p matrix of covariates, and the rows are i.i.d multivariate normal N(0, Omega^-1).
We find sparse estimates of B and Omega using a multivariate extension of the spike-and-slab LASSO (SSL) regularization found in .
Specifically, we place independent SSL priors on the entries of the matrix of covariate effects B and the off-diagonal elements of the residual precision matrix Omega. 
\code{mSSL_dpe} runs an Expectation-Conditional Maximization algorithm that targets the MAP estimates of B and Omega over a grid of posteriors 
corresponding to a user-specified sequence of spike distributions for both the covariate effects and partial covariances. 
See Deshpande et al. (2018) for details of how the function navigates the grid of posteriors.
Note: between \code{\link{mSSL_dcpe}} and \code{\link{sep_sslg}}, \code{mSSL_dpe} is the slowest but the most accurate.
We recommend running both \code{mSSL_dpe} and \code{\link{mSSL_dcpe}}.}

\usage{
mSSL_dpe(X, Y, 
         lambdas = list(lambda1 = 1, lambda0 = seq(1, nrow(X), length = 10)), 
         xis = list(xi1 = 0.01 * nrow(X),  xi0 = seq(0.1 * nrow(X), nrow(X), length = 10)), 
         theta_hyper_params = c(1, ncol(X) * ncol(Y)), 
         eta_hyper_params = c(1, ncol(Y)), 
         diag_penalty = 1,
         max_iter = 500,
         eps = 1e-3,
         s_max_condition = 10*nrow(X),
         obj_counter_max = 5,
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
  \item{diag_penalty}{Binary indicator whether or not the diagonal elements of Omega should be penalized (\code{diag_penalty = 1}). Default is 1.}
  \item{max_iter}{Maximum number of iterations for the main ECM algorithm. Default value is 500.}
  \item{eps}{Convergence tolerance for the ECM algorithm. If percetage change in every element of B and Omega is less than \code{eps}, the algorithm terminates. Default value is 0.001.}
  \item{s_max_condition}{Maximum allowable condition number of the residual covariance matrix S = (Y - XB)'(Y - XB). If condition number of S exceeds \code{s_max_condition}, the algorithm terminates. This typically occurs only if B is dense and there is very little residual variation. Default value is 10 * n}
  \item{obj_counter_max}{When the percentage increase in objective function is less than \code{eps} for \code{obj_counter_max} consecutive iterations, algorithm terminates. Default value is 5.}
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
  \item{B_path}{A p x q x L^2 array containing the MAP estimate of B corresponding to each posterior explored}
  \item{Omega_path}{A q x q x L^2 array containing the MAP estimate of Omega corresponding to each posterior explored}
  \item{theta_path}{An L x L matrix containing the MAP estimate of theta from each posterior explored}
  \item{eta_path}{An L x L matrix containing the MAP estimate of eta from each posterior explored}
  \item{time}{Total time in seconds}
  \item{lambda0}{The vector of lambda0's used. This is used for plotting purposes. See \code{\link{plot_support}}}
  \item{xi0}{The vector of xi0's used. This is primarily used for plotting purposes. See \code{\link{plot_support}}}
  \item{early_term}{A vector of length L^2 of binary indicators that tell whether the ECM algorithm terminated early due to a severely ill-conditioned residual covariance matrix. This is primarily used for plotting purposes. See \code{\link{plot_support}}.}
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
fit_mSSL_dpe <- mSSL_dpe(X, Y)

}
