\name{error_B}
\alias{error_B}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Function to compute various error measures for B
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function computes various performance measures to evaluate the estimated B. It includes measures of support recovery as well as estimation accuracy. Note: using this function only makes sense for simulation experiments where we know the actual B that generated the data
}
\usage{
error_B(B, B_orig)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{B}{Estimated matrix B
%%     ~~Describe \code{x} here~~
}
  \item{B_orig}{The B matrix used to generate the data}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x) 
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
