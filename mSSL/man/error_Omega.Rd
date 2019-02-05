\name{error_Omega}
\alias{error_Omega}
\title{
Function to compute various error measures for B
}
\description{
This function computes various performance measures to evaluate the estimated Omega. It includes measures of support recovery as well as estimation accuracy. Note: using this function only makes sense for simulation experiments where we know the actual Omega that generated the data
}
\usage{
error_Omega(Omega, Omega_orig)
}
\arguments{
  \item{Omega}{Estimated matrix Omega}
  \item{Omega_orig}{The Omega matrix used to generate the data}
}
\details{
}
\value{
Returns a vector with the following named entries
  \item{TP}{Number of true positives}
  \item{TN}{Number of true negatives}
  \item{FP}{Number of false positives}
  \item{FN}{Number of false negatives}
  \item{SEN}{Sensitivity. Equal to TP/(TP + FN)}
  \item{SPE}{Specificity. Equal to TN/(TN + FP)}
  \item{PREC}{Precision. Equal to TP/(TP + TN)}
  \item{ACC}{Accuracy. Equal to (TP + TN)/(TP + TN + FP + FN)}
  \item{F1}{F1 Score. Equal to (2 * TP)/(2 * TP + FP + FN)}
  \item{MCC}{Matthew Correlation Coefficient. Equal to (TP/N - S * P)/sqrt(PS(1 - S)(1 - P)), where N = TN + TP + FN + FP, S = (TP + FN)/N, and P = (TP + FP)/N.}
  \item{FROB}{Frobenius norm of Omega - Omega_orig.}
}
\references{
}
\author{
Sameer Deshpande
}
\note{
}


\seealso{
\code{\link{error_B}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
