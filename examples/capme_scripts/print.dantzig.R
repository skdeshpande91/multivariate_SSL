print.dantzig <- function(x,digits = max(3, getOption("digits") - 3), ... ) {
  cat("\n lambdas used:\n")
  print(signif(x$lambda,digits))
}
