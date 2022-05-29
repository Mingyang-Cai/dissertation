#SMC FCS interaction term
rm(list=ls())
library(mice, warn.conflicts = FALSE)
library(miceadds)
library(MASS)
library(smcfcs)

set.seed(123)

sample.size <- 10000
missingness <- 0.5 #percentage Y incomplete
simulations <- 100

#prepare output object
coefficient.OUT <- list()

#start simulation
pb <- txtProgressBar(min = 0, max = simulations, style = 3)
for (i in 1:simulations){
  # Generate data
  x1 <- rnorm(sample.size, mean = 0, sd = 1)
  x2 <- rnorm(sample.size, mean = 0, sd = 1)
  epsilon <- rnorm(sample.size, mean = 0, sd = 1)
  y <- x1 + x2 + x1*x2 + epsilon
  
  complete.data <- cbind(y, x1, x2,  x1*x2)
  
  
  # Generate missings in Y
  incomplete.data <- ampute(data = complete.data, prop = missingness, mech = "MCAR", patterns = c(1, 0, 0, 0))$amp
  
  #incomplete.data <- ampute(data=complete.data, prop=missingness, mech = "MAR", type = 'LEFT', patterns = c(1, 0, 0, 0),
  #                          weights = c(1, 0, 0, 0))$amp
  colnames(incomplete.data) <- c("Y", "X1", "X2", "X1X2")
  
  
  imp <- smcfcs(incomplete.data, smtype = "lm", smformula = "Y~X1+X2+X1X2", method = c("", "norm", "norm", "X1*X2"), rjlimit= 1000)
  
  summary.OUT <- list()
  for (j in 1:5) {
    summary.OUT[[j]] <- summary(lm(imp$impDatasets[[j]]$Y~imp$impDatasets[[j]]$X1+imp$impDatasets[[j]]$X2+imp$impDatasets[[j]]$X1X2))
  }
  
  coefficient.OUT[[i]] <- summary.OUT
  
  
  setTxtProgressBar(pb, i)
}
close(pb)


# Evaluate the coverage rate
evaluate.sims <- function(sims, truth = 1){
  POOL <- list()
  pb <- txtProgressBar(min = 0, max = simulations, style = 3)
  for (i in 1:length(sims)){
    #Extract means and variances
    Q            <- c(coefficient.OUT[[i]][[1]]$coefficients[4,1], coefficient.OUT[[i]][[2]]$coefficients[4,1],
                      coefficient.OUT[[i]][[3]]$coefficients[4,1], coefficient.OUT[[i]][[4]]$coefficients[4,1],
                      coefficient.OUT[[i]][[5]]$coefficients[4,1])
    
    U            <- c((coefficient.OUT[[i]][[1]]$coefficients[4,2])^2, (coefficient.OUT[[i]][[2]]$coefficients[4,2])^2,
                      (coefficient.OUT[[i]][[3]]$coefficients[4,2])^2, (coefficient.OUT[[i]][[4]]$coefficients[4,2])^2, 
                      (coefficient.OUT[[i]][[5]]$coefficients[4,2])^2)
    
    #Pool the regular way
    pool         <- mice::pool.scalar(Q, U, n = 1000) # A really large number
    pool$lower   <- pool$qbar - qt(0.975, pool$df) * sqrt(pool$t)
    pool$upper   <- pool$qbar + qt(0.975, pool$df) * sqrt(pool$t)
    pool$coverage <- pool$lower < mean(truth) & mean(truth) < pool$upper
    POOL[[i]]     <- unlist(pool)
    setTxtProgressBar(pb, i)
  }
  return(POOL)
  close(pb)
}
EVAL <- evaluate.sims(coefficient.OUT)

# Summarize
AVG.EVAL <- Reduce("+", EVAL) / length(EVAL)
round(AVG.EVAL, 3)

# Evaluate the coverage rate
evaluate.sims <- function(sims, truth = 1){
  POOL <- list()
  pb <- txtProgressBar(min = 0, max = simulations, style = 3)
  for (i in 1:length(sims)){
    #Extract means and variances
    Q            <- c(coefficient.OUT[[i]][[1]]$coefficients[3,1], coefficient.OUT[[i]][[2]]$coefficients[3,1],
                      coefficient.OUT[[i]][[3]]$coefficients[3,1], coefficient.OUT[[i]][[4]]$coefficients[3,1],
                      coefficient.OUT[[i]][[5]]$coefficients[3,1])
    
    U            <- c((coefficient.OUT[[i]][[1]]$coefficients[3,2])^2, (coefficient.OUT[[i]][[2]]$coefficients[3,2])^2,
                      (coefficient.OUT[[i]][[3]]$coefficients[3,2])^2, (coefficient.OUT[[i]][[4]]$coefficients[3,2])^2, 
                      (coefficient.OUT[[i]][[5]]$coefficients[3,2])^2)
    
    #Pool the regular way
    pool         <- mice::pool.scalar(Q, U, n = 1000) # A really large number
    pool$lower   <- pool$qbar - qt(0.975, pool$df) * sqrt(pool$t)
    pool$upper   <- pool$qbar + qt(0.975, pool$df) * sqrt(pool$t)
    pool$coverage <- pool$lower < mean(truth) & mean(truth) < pool$upper
    POOL[[i]]     <- unlist(pool)
    setTxtProgressBar(pb, i)
  }
  return(POOL)
  close(pb)
}
EVAL <- evaluate.sims(coefficient.OUT)

# Summarize
AVG.EVAL <- Reduce("+", EVAL) / length(EVAL)
round(AVG.EVAL, 3)






# Evaluate the coverage rate
evaluate.sims <- function(sims, truth = 1){
  POOL <- list()
  pb <- txtProgressBar(min = 0, max = simulations, style = 3)
  for (i in 1:length(sims)){
    #Extract means and variances
    Q            <- c(coefficient.OUT[[i]][[1]]$coefficients[2,1], coefficient.OUT[[i]][[2]]$coefficients[2,1],
                      coefficient.OUT[[i]][[3]]$coefficients[2,1], coefficient.OUT[[i]][[4]]$coefficients[2,1],
                      coefficient.OUT[[i]][[5]]$coefficients[2,1])
    
    U            <- c((coefficient.OUT[[i]][[1]]$coefficients[2,2])^2, (coefficient.OUT[[i]][[2]]$coefficients[2,2])^2,
                      (coefficient.OUT[[i]][[3]]$coefficients[2,2])^2, (coefficient.OUT[[i]][[4]]$coefficients[2,2])^2, 
                      (coefficient.OUT[[i]][[5]]$coefficients[2,2])^2)
    
    #Pool the regular way
    pool         <- mice::pool.scalar(Q, U, n = 1000) # A really large number
    pool$lower   <- pool$qbar - qt(0.975, pool$df) * sqrt(pool$t)
    pool$upper   <- pool$qbar + qt(0.975, pool$df) * sqrt(pool$t)
    pool$coverage <- pool$lower < mean(truth) & mean(truth) < pool$upper
    POOL[[i]]     <- unlist(pool)
    setTxtProgressBar(pb, i)
  }
  return(POOL)
  close(pb)
}
EVAL <- evaluate.sims(coefficient.OUT)

# Summarize
AVG.EVAL <- Reduce("+", EVAL) / length(EVAL)
round(AVG.EVAL, 3)



intercept<-rep(0, 100)
for (i in 1:100) {
  intercept[i] <- mean(coefficient.OUT[[i]][[1]]$coefficients[1,1], coefficient.OUT[[i]][[2]]$coefficients[1,1], coefficient.OUT[[i]][[3]]$coefficients[1,1],
                       coefficient.OUT[[i]][[4]]$coefficients[1,1], coefficient.OUT[[i]][[5]]$coefficients[1,1])
}

mean(intercept)


resi<-rep(0, 100)
for (i in 1:100) {
  resi[i] <- mean(coefficient.OUT[[i]][[1]]$sigma, coefficient.OUT[[i]][[2]]$sigma, coefficient.OUT[[i]][[3]]$sigma,
                  coefficient.OUT[[i]][[4]]$sigma, coefficient.OUT[[i]][[5]]$sigma)
}

mean(resi)


rsquare<-rep(0, 100)
for (i in 1:100) {
  rsquare[i] <- mean(coefficient.OUT[[i]][[1]]$r.square, coefficient.OUT[[i]][[2]]$r.square, coefficient.OUT[[i]][[3]]$r.square,
                     coefficient.OUT[[i]][[4]]$r.square, coefficient.OUT[[i]][[5]]$r.square)
}

mean(rsquare)