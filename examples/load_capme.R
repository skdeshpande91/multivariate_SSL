# Load all of the CAPME scripts
source("scripts/capme_scripts/capme.R")
source("scripts/capme_scripts/cv.capme.R")
source("scripts/capme_scripts/cv.part.R")
source("scripts/capme_scripts/dantzig.R")
source("scripts/capme_scripts/likelihood.R")
source("scripts/capme_scripts/linprogPD2.R")
source("scripts/capme_scripts/linprogS2.R")
source("scripts/capme_scripts/make.pd.R")
source("scripts/capme_scripts/meanclime.R")
source("scripts/capme_scripts/Mpower.R")
source("scripts/capme_scripts/print.capme.R")
source("scripts/capme_scripts/print.cv.capme.R")
source("scripts/capme_scripts/print.dantzig.R")
source("scripts/capme_scripts/tracel2.R")
source("scripts/capme_scripts/zzz.R")

my.capme <- function(X, Y, linsolver.Gamma = c("simplex", "primaldual"), linsolver.Omega = c("simplex", "primaldual")){
  tmp.fit <- cv.capme(fold = 5, loss = "likelihood", X, Y, linsolver.Gamma = linsolver.Gamma, linsolver.Omega = linsolver.Omega)
  lambdaopt <- tmp.fit$lambdaopt
  tauopt <- tmp.fit$tauopt
  capme.obj <- capme(X, Y, lambda = lambdaopt, tau = tauopt, linsolver.Gamma = linsolver.Gamma, linsolver.Omega = linsolver.Omega)
  results <- list("B" = capme.obj$Gammalist[[1]], "Omega" = capme.obj$Omegalist[[1]])
  return(results)
}

# trying capme

#test.capme <- cv.capme(fold = 5, loss = "likelihood", X, Y, linsolver.Gamma = "simplex", linsolver.Omega = "simplex")
# Now we need to resolve the object
#lambdaopt <- test.capme$lambdaopt
#tauopt <- test.capme$tauopt

#capme.obj <- capme(X, Y, lambda = lambdaopt, tau = tauopt, linsolver.Gamma = "simplex", linsolver.Omega = "simplex")
#capme.obj$Gammalist
#capme.obj$Omegalist
