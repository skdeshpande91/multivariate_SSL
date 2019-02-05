error_B <- function(B, B_orig){
  
    perf <- rep(NA, times = 11)
    names(perf) <- c("TP", "TN", "FP", "FN", "SEN", "SPE", "PREC", "ACC",  "F1", "MCC", "MSE")
    tp <- sum(B != 0 & B_orig != 0)
    tn <- sum(B == 0 & B_orig == 0)
    fp <- sum(B != 0 & B_orig == 0)
    fn <- sum(B == 0 & B_orig != 0)
    
    sen <- tp/(tp + fn)
    spe <- tn/(tn + fp)
    prec <- tp/(tp + fp)
    acc <- (tp + tn)/(tp + tn + fp + fn)
    f1 <- 2 * tp/(2*tp + fp + fn)
    
    mcc.n <- tp + tn + fp + fn
    mcc.s <- (tp + fn)/mcc.n
    mcc.p <- (tp + fp)/mcc.n
    mcc <- (tp/mcc.n - mcc.s*mcc.p)/sqrt(mcc.p * mcc.s * (1 - mcc.p) * (1 - mcc.s))
    perf["TP"] <- tp 
    perf["TN"] <- tn 
    perf["FP"] <- fp 
    perf["FN"] <- fn 
    perf["SEN"] <- sen
    perf["SPE"] <- spe
    perf["PREC"] <- prec
    perf["ACC"] <- acc
    perf["F1"] <- f1
    perf["MCC"] <- mcc
    perf["MSE"] <- mean((B - B_orig)^2)
  return(perf)

}
