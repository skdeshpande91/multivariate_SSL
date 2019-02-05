print.capme <- function(x,digits = max(3, getOption("digits") - 3), ... )
{
  cat(" perturb=", signif(x$perturb, digits))
  cat("\n lambds used: \n")
  print(signif(x$lambda,digits))
  cat("\n taus used:\n")
  print(signif(x$tau,digits))
}
