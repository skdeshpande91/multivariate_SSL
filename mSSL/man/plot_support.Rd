\name{plot_support}
\alias{plot_support}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Stability Plot of the Estimated Supports of B and Omega
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Visual self-check of whether the supports of the estimated B and Omega have stabilized along the mSSL-DPE trajectory.
}
\usage{
plot_support(fit_mSSL_dpe)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit_mSSL_dpe}{Object returned by a call to \code{\link{mSSL_dpe}}}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
There is a point in the plot for each combination of lambda0 and xi0 along the DPE path.
Red X's indicate posterior means which may be numerical unstable (i.e. have residual covariance matrices with very high condition number). 
Blue +'s indicate that the support of the estimated B and Omega for the corresponding combination of spike penalties is identical to the estimated supports at the end of the path.
Black dots correspond to modes where the estimated supports differ from the final supports.
}
\value{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
