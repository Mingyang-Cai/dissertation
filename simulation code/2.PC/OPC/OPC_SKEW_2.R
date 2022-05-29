#OPC
rm(list=ls())
library(mice, warn.conflicts = FALSE)
library(miceadds)
library(MASS)
library(sn)

set.seed(123)

#prepare output object
coefficient.OUT <- list()

Ampute.function <- function(data, missingness, m.mech = "MARleft"){
  incomplete.data <- data
  if(m.mech == "MARleft"){
    incomplete.data <- ampute(data = data, prop = missingness, mech = "MAR", type = 'LEFT', patterns = c(0, 0, 1),
                              weights = c(0, 0, 1))$amp
  }
  
  if(m.mech == "MARmid"){
    incomplete.data <- ampute(data = data, prop = missingness, mech = "MAR", type = 'MID', patterns = c(0, 0, 1),
                              weights = c(0, 0, 1))$amp
  }
  
  if(m.mech == "MARtail"){
    incomplete.data <- ampute(data = data, prop = missingness, mech = "MAR", type = 'TAIL', patterns = c(0, 0, 1),
                              weights = c(0, 0, 1))$amp
  }
  
  if(m.mech == "MARright"){
    incomplete.data <- ampute(data = data, prop=missingness, mech = "MAR", type = 'RIGHT', patterns = c(0, 0, 1),
                              weights = c(0, 0, 1))$amp
  }
  
  if(m.mech == "MCAR"){
    incomplete.data <- ampute(data = data, prop = missingness, mech = "MCAR", patterns = c(0, 0, 1))$amp
  }
  return(incomplete.data)
}

OPC.function <- function(complete.data, missingness, m.mech){
  incomplete.data <- Ampute.function(complete.data, missingness, m.mech)
  ini <- mice(incomplete.data, maxit = 0)
  meth <- c("quadratic", "~I(X^2)", "")
  pred <- ini$pred
  pred[, "XX"] <- 0
  imp.opc <- mice(incomplete.data, meth = meth, pred = pred, print = FALSE)
  fit <- with(imp.opc, lm(Y ~ X + XX))
  pool <- pool(fit)$pooled
  pool$lower <- pool$estimate - (qt(.975, pool$df) * sqrt(pool$t))
  pool$upper <- pool$estimate + (qt(.975, pool$df) * sqrt(pool$t))
  pool$true <- c(0, 1, 1)
  pool$bias <- pool$estimate - pool$true
  pool$cov <- pool$lower < pool$true & pool$true < pool$upper
  pool$ciw <- pool$upper - pool$lower
  return(pool)
}

sample.size <- 100
missingness <- 0.3 #percentage X incomplete
simulations <- 1000



MCAR <- list()
MARl <- list()
MARm <- list()
MARt <- list()
MARr <- list()

pb <- txtProgressBar(min = 0, max = simulations, style = 3)
for (i in 1 : simulations) {
  X <- rsn(sample.size, xi = 2, omega = 1, alpha = -3)
  epsilon <- rnorm(sample.size, mean = 0, sd = 1.2)
  Y <- X + X^2 + epsilon
  complete.data <- data.frame(
    X,
    XX = X^2,
    Y
  )
  MCAR[[i]] <-OPC.function(complete.data = complete.data, missingness = missingness, m.mech = "MCAR")
  MARl[[i]] <-OPC.function(complete.data = complete.data, missingness = missingness, m.mech = "MARleft")
  MARm[[i]] <-OPC.function(complete.data = complete.data, missingness = missingness, m.mech = "MARmid")
  MARt[[i]] <-OPC.function(complete.data = complete.data, missingness = missingness, m.mech = "MARtail")
  MARr[[i]] <-OPC.function(complete.data = complete.data, missingness = missingness, m.mech = "MARright")
  setTxtProgressBar(pb, i)
}
close(pb) 
res.MCAR <- round(Reduce("+", MCAR) / length(MCAR), 2)
res.MARl <- round(Reduce("+", MARl) / length(MARl), 2)
res.MARm <- round(Reduce("+", MARm) / length(MARm), 2)
res.MARt <- round(Reduce("+", MARt) / length(MARt), 2)
res.MARr <- round(Reduce("+", MARr) / length(MARr), 2)

