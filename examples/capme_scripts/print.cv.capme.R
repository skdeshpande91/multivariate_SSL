print.cv.capme <- function(x,digits = max(3, getOption("digits") - 3), ... )
{
  cat("\n lambds used: \n")
  print(signif(x$lambda,digits))
  cat("\n taus used:\n")
  print(signif(x$tau,digits))
  cat("\n average loss: (rows: lambda; columns: tau) \n")
  print(signif(x$loss.mean,digits))
  cat("\n lambdaopt: \n")
  print(signif(x$lambdaopt,digits))
  cat("\n tauopt: \n")
  print(signif(x$tauopt,digits))
}
