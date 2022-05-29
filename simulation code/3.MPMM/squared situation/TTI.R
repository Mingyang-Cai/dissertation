## TTI  COVER

rm(list=ls())
library(mice, warn.conflicts = FALSE)
library(miceadds)
library(MASS)

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
  y <- rnorm(sample.size, mean = 0, sd = 1)
  epsilon <- rnorm(sample.size, mean = 0, sd = 1)
  x <- y + y^2 + epsilon
  
  complete.data <- cbind(y, y^2, x)
  
  
  # Generate missings in Y
  #incomplete.data <- ampute(data = complete.data, prop = missingness, mech = "MCAR", patterns = c(0, 0, 1))$amp
  
  incomplete.data <- ampute(data=complete.data, prop=missingness, mech = "MAR", type = 'RIGHT', patterns = c(0, 0, 1),
                            weights = c(0, 0, 1))$amp
  colnames(incomplete.data) <- c("Y", "YY", "X")
  
  data <- incomplete.data
  imp <- mice(data, m = 5, method = "norm",  print = FALSE)
  
  summary.OUT <- list()
  for (j in 1:5) {
    summary.OUT[[j]] <- summary(lm(complete(imp,j)$X~complete(imp,j)$Y+complete(imp,j)$YY))
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

