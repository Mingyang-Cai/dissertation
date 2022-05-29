## PMM with CRA

rm(list=ls())
library(mice, warn.conflicts = FALSE)
library(miceadds)
library(MASS)

set.seed(123)

sample.size <- 100
missingness <- 0.3 #percentage Y incomplete
simulations <- 1000






pmm.cra.single <- function(data){
  data.mis <- data[is.na(data[, 1])==TRUE, ]
  data.obs <- data[is.na(data[, 1])==FALSE, ]
  
  data.obs.y <- data.obs[1:2]
  data.obs.x <- data.obs[3]
  
  data.mis.y <- data.mis[1:2]
  data.mis.x <- data.mis[3]
  
  ES <- eigen(solve(cov(data.obs.y, data.obs.y)) %*% cov(data.obs.y, data.obs.x) 
              %*% solve(cov(data.obs.x, data.obs.x)) %*% cov(data.obs.x, data.obs.y))
  
  alpha <- ES$vectors[,1]
  alpha <- as.matrix(alpha)
  
  
  frv.obs <- as.matrix(data.obs.y) %*% alpha
  frv.mis <- as.matrix(data.mis.y) %*% alpha
  
  
  
  frv.total <- rbind(frv.obs, frv.mis)
  #test.data <- cbind(frv.total, rbind(data.obs.x, data.mis.x))
  
  imp<-mice.impute.pmm(frv.total,!is.na(frv.total),rbind(data.obs.x, data.mis.x), donors = 10L)
  
  frv.mis.est <- as.matrix(imp)
  
  for (i in 1:dim(frv.mis.est)[1]){
    
    data.mis.y[i,]<-data.obs.y[which(frv.obs==frv.mis.est[i])[1],]
  }
  
  cbind(rbind(data.obs.y, data.mis.y), rbind(data.obs.x, data.mis.x))
}

pmm.cra.multiple <- function(data){
  multiple.imputed.dataset <- list()
  
  for(i in 1:5){
    multiple.imputed.dataset[[i]] <- pmm.cra.single(data)
  }
  
  multiple.imputed.dataset
}




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
  incomplete.data <- ampute(data = complete.data, prop = missingness, mech = "MCAR", patterns = c(0, 0, 1))$amp
  
  #incomplete.data <- ampute(data=complete.data, prop=missingness, mech = "MAR", type = 'RIGHT', patterns = c(0, 0, 1),
  #                         weights = c(0, 0, 1))$amp
  colnames(incomplete.data) <- c("Y", "YY", "X")
  
  
  
  
  # Impute
  OUT <- pmm.cra.multiple(incomplete.data)
  
  
  summary.OUT <- list()
  for (j in 1:5) {
    summary.OUT[[j]] <- summary(lm(OUT[[j]]$X~OUT[[j]]$Y+OUT[[j]]$YY))
  }
  
  coefficient.OUT[[i]]<- summary.OUT
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

