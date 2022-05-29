#NIMIP
rm(list=ls())
library(mice, warn.conflicts = FALSE)
library(miceadds)
library(MASS)
library(magrittr)
library(purrr)

nsim <- 1000
rho.m <- 0.8                               #marginal correlation
#rho.p <- (0.8 -0.5*0.5) / (1 - 0.5^2)      #partial correlation
#rho.p <- 0
sample.size <- 5000
set.seed(123)
coefficient.mean <- list()
coefficient.var <- list()
#fisher's transformation
fisher.trans <- function(x){
  z <- log((1 + x) / (1 - x)) / 2
  return(z)
}
#create a single imputed dataset
impute.single <- function(incomplete.data, rho.p){
parm.y0 <- as.matrix(incomplete.data$x) %>% 
          cbind(1, .) %>% 
          .norm.draw(incomplete.data$y0, !is.na(incomplete.data$y0), .)
parm.y1 <- as.matrix(incomplete.data$x) %>% 
          cbind(1, .) %>% 
          .norm.draw(incomplete.data$y1, !is.na(incomplete.data$y1), .)
cov.y0y1.p <- rho.p * parm.y0$sigma * parm.y1$sigma
#impute y0
x0 <- as.matrix(incomplete.data$x) %>% 
  cbind(1, .) %>% 
  extract(is.na(incomplete.data$y0), )
imp.y0 <- x0 %>% 
          multiply_by_matrix(parm.y0$beta) %>%
          add((as.matrix(incomplete.data$y1[is.na(incomplete.data$y0)]) - x0 %*% parm.y1$beta) * cov.y0y1.p / (parm.y1$sigma)^2) %>%
          add(rnorm(sum(is.na(incomplete.data$y0))) * sqrt((1 - rho.p^2)) * parm.y0$sigma)
#impute y1
x1 <- as.matrix(incomplete.data$x) %>% 
  cbind(1, .) %>% 
  extract(is.na(incomplete.data$y1), )
imp.y1 <- x1 %>% 
  multiply_by_matrix(parm.y1$beta) %>%
  add((as.matrix(incomplete.data$y0[is.na(incomplete.data$y1)]) - x1 %*% parm.y0$beta) * cov.y0y1.p / (parm.y0$sigma)^2) %>%
  add(rnorm(sum(is.na(incomplete.data$y1))) * sqrt((1 - rho.p^2)) * parm.y1$sigma)
incomplete.data$y0[1 : 2500] <- imp.y0
incomplete.data$y1[2501 : 5000] <- imp.y1
return(incomplete.data)
}
#create multiple imputed datasets
impute.multiple <- function(m = 5, incomplete.data, rho.p, ...){
 multi.imp.dat <- list()
 for (i in 1 : m) {
 multi.imp.dat[[i]] <-impute.single(incomplete.data, rho.p)  
 }
 return(multi.imp.dat)
}
pb <- txtProgressBar(min = 0, max = nsim, style = 3)
for (i in 1 : nsim) {
  data <- mvrnorm(n = sample.size, mu = c(0, 1, 2), 
                  Sigma = matrix(c(1, rho.m, 0.5, rho.m, 1, 0.5, 0.5, 0.5, 1), nrow = 3)) %>% as.data.frame
  colnames(data)<- c("y0", "y1", "x")
  incomplete.data <- data 
  incomplete.data$y0[1 : 2500] <- NA
  incomplete.data$y1[2501 : 5000] <- NA
  OUT <- impute.multiple(5, incomplete.data, (0.8 -0.5*0.5) / (1 - 0.5^2))
  mean.OUT <- list()
  variance.OUT <- list()
  for (j in 1 : 5) {
    mean.OUT[[j]] <- colMeans(OUT[[j]])
    variance.OUT[[j]] <- cov(OUT[[j]])
  }
  coefficient.mean[[i]] <- mean.OUT
  coefficient.var[[i]] <- variance.OUT
  setTxtProgressBar(pb, i)
}
close(pb)


# Evaluate mean of y0
truth <- 0  
POOL <- list()
  for (i in 1:nsim){
    #Extract means and variances
    Q            <- c(coefficient.mean[[i]][[1]][1], coefficient.mean[[i]][[2]][1],
                      coefficient.mean[[i]][[3]][1], coefficient.mean[[i]][[4]][1],
                      coefficient.mean[[i]][[5]][1])
    
    U            <- c(coefficient.var[[i]][[1]][1,1], coefficient.var[[i]][[2]][1,1],
                      coefficient.var[[i]][[3]][1,1], coefficient.var[[i]][[4]][1,1], 
                      coefficient.var[[i]][[5]][1,1]) / sample.size
    #Pool the regular way
    pool         <- mice::pool.scalar(Q, U, n = 1000) # A really large number
    pool$lower   <- pool$qbar - qt(0.975, pool$df) * sqrt(pool$t)
    pool$upper   <- pool$qbar + qt(0.975, pool$df) * sqrt(pool$t)
    pool$coverage <- pool$lower < mean(truth) & mean(truth) < pool$upper
    POOL[[i]]     <- unlist(pool)
  }
AVG.mean.y0 <- Reduce("+", POOL) / length(POOL)
round(AVG.mean.y0, 3)

# Evaluate mean of y1
truth <- 1  
POOL <- list()
for (i in 1:nsim){
  #Extract means and variances
  Q            <- c(coefficient.mean[[i]][[1]][2], coefficient.mean[[i]][[2]][2],
                    coefficient.mean[[i]][[3]][2], coefficient.mean[[i]][[4]][2],
                    coefficient.mean[[i]][[5]][2])
  
  U            <- c(coefficient.var[[i]][[1]][2,2], coefficient.var[[i]][[2]][2,2],
                    coefficient.var[[i]][[3]][2,2], coefficient.var[[i]][[4]][2,2], 
                    coefficient.var[[i]][[5]][2,2]) / sample.size
  #Pool the regular way
  pool         <- mice::pool.scalar(Q, U, n = 1000) # A really large number
  pool$lower   <- pool$qbar - qt(0.975, pool$df) * sqrt(pool$t)
  pool$upper   <- pool$qbar + qt(0.975, pool$df) * sqrt(pool$t)
  pool$coverage <- pool$lower < mean(truth) & mean(truth) < pool$upper
  POOL[[i]]     <- unlist(pool)
}
AVG.mean.y1 <- Reduce("+", POOL) / length(POOL)
round(AVG.mean.y1, 3)

# Evaluate var of y0
truth <- 1  
POOL <- list()
for (i in 1:nsim){
  #Extract means and variances
  Q            <- c(coefficient.var[[i]][[1]][1,1], coefficient.var[[i]][[2]][1,1],
                    coefficient.var[[i]][[3]][1,1], coefficient.var[[i]][[4]][1,1],
                    coefficient.var[[i]][[5]][1,1])
  
  U            <- c(coefficient.var[[i]][[1]][1,1]^2, coefficient.var[[i]][[2]][1,1]^2,
                    coefficient.var[[i]][[3]][1,1]^2, coefficient.var[[i]][[4]][1,1]^2,
                    coefficient.var[[i]][[5]][1,1]^2) * 2 / sample.size
  #Pool the regular way
  pool         <- mice::pool.scalar(Q, U, n = 1000) # A really large number
  pool$lower   <- pool$qbar - qt(0.975, pool$df) * sqrt(pool$t)
  pool$upper   <- pool$qbar + qt(0.975, pool$df) * sqrt(pool$t)
  pool$coverage <- pool$lower < mean(truth) & mean(truth) < pool$upper
  POOL[[i]]     <- unlist(pool)
}
AVG.var.y0 <- Reduce("+", POOL) / length(POOL)
round(AVG.var.y0, 3)

# Evaluate var of y1
truth <- 1  
POOL <- list()
for (i in 1:nsim){
  #Extract means and variances
  Q            <- c(coefficient.var[[i]][[1]][2,2], coefficient.var[[i]][[2]][2,2],
                    coefficient.var[[i]][[3]][2,2], coefficient.var[[i]][[4]][2,2],
                    coefficient.var[[i]][[5]][2,2])
  
  U            <- c(coefficient.var[[i]][[1]][2,2], coefficient.var[[i]][[2]][2,2],
                    coefficient.var[[i]][[3]][2,2], coefficient.var[[i]][[4]][2,2],
                    coefficient.var[[i]][[5]][2,2]) * 2 / sample.size
  #Pool the regular way
  pool         <- mice::pool.scalar(Q, U, n = 1000) # A really large number
  pool$lower   <- pool$qbar - qt(0.975, pool$df) * sqrt(pool$t)
  pool$upper   <- pool$qbar + qt(0.975, pool$df) * sqrt(pool$t)
  pool$coverage <- pool$lower < mean(truth) & mean(truth) < pool$upper
  POOL[[i]]     <- unlist(pool)
}
AVG.var.y1 <- Reduce("+", POOL) / length(POOL)
round(AVG.var.y1, 3)

# Evaluate cov of y0y1
truth <- 0.8  
POOL <- list()
for (i in 1:nsim){
  #Extract means and variances
  Q            <- c(coefficient.var[[i]][[1]][1,2]/sqrt(coefficient.var[[i]][[1]][1,1] * coefficient.var[[i]][[1]][2,2]), 
                    coefficient.var[[i]][[2]][1,2]/sqrt(coefficient.var[[i]][[2]][1,1] * coefficient.var[[i]][[2]][2,2]),
                    coefficient.var[[i]][[3]][1,2]/sqrt(coefficient.var[[i]][[3]][1,1] * coefficient.var[[i]][[3]][2,2]), 
                    coefficient.var[[i]][[4]][1,2]/sqrt(coefficient.var[[i]][[4]][1,1] * coefficient.var[[i]][[4]][2,2]),
                    coefficient.var[[i]][[5]][1,2]/sqrt(coefficient.var[[i]][[5]][1,1] * coefficient.var[[i]][[5]][2,2])) %>% 
                    fisher.trans()
  U            <- rep(1 / (sample.size - 3), 5)
  #Pool the regular way
  pool         <- mice::pool.scalar(Q, U, n = 1000) # A really large number
  lower.bound  <- pool$qbar - qt(0.975, pool$df) * sqrt(pool$t)
  upper.bound  <- pool$qbar + qt(0.975, pool$df) * sqrt(pool$t)
  
  
  pool$lower   <- (exp(2 * lower.bound) - 1) / (exp(2 * lower.bound) + 1)
  pool$upper   <- (exp(2 * upper.bound) - 1) / (exp(2 * upper.bound) + 1)
  pool$coverage <- pool$lower < mean(truth) & mean(truth) < pool$upper
  pool$qbar <- (exp(2 * pool$qbar) - 1) / (exp(2 * pool$qbar) + 1)
  POOL[[i]]     <- c(pool$m, pool$qbar, pool$lower, pool$upper, pool$coverage)
}
AVG.cov.y0y1 <- Reduce("+", POOL) / length(POOL)
names(AVG.cov.y0y1) <- c("m", "qbar", "lower", "upper", "coverage")
round(AVG.cov.y0y1, 3)

# Evaluate cov of y0x
truth <- 0.5  
POOL <- list()
for (i in 1:nsim){
  #Extract means and variances
  Q            <- c(coefficient.var[[i]][[1]][1,3]/sqrt(coefficient.var[[i]][[1]][1,1] * coefficient.var[[i]][[1]][3,3]), 
                    coefficient.var[[i]][[2]][1,3]/sqrt(coefficient.var[[i]][[2]][1,1] * coefficient.var[[i]][[2]][3,3]),
                    coefficient.var[[i]][[3]][1,3]/sqrt(coefficient.var[[i]][[3]][1,1] * coefficient.var[[i]][[3]][3,3]), 
                    coefficient.var[[i]][[4]][1,3]/sqrt(coefficient.var[[i]][[4]][1,1] * coefficient.var[[i]][[4]][3,3]),
                    coefficient.var[[i]][[5]][1,3]/sqrt(coefficient.var[[i]][[5]][1,1] * coefficient.var[[i]][[5]][3,3])) %>% 
                    fisher.trans()
  
  U            <- rep(1 / (sample.size - 3), 5)
  #Pool the regular way
  pool         <- mice::pool.scalar(Q, U, n = 1000) # A really large number
  lower.bound  <- pool$qbar - qt(0.975, pool$df) * sqrt(pool$t)
  upper.bound  <- pool$qbar + qt(0.975, pool$df) * sqrt(pool$t)
  
  
  pool$lower   <- (exp(2 * lower.bound) - 1) / (exp(2 * lower.bound) + 1)
  pool$upper   <- (exp(2 * upper.bound) - 1) / (exp(2 * upper.bound) + 1)
  pool$coverage <- pool$lower < mean(truth) & mean(truth) < pool$upper
  pool$qbar <- (exp(2 * pool$qbar) - 1) / (exp(2 * pool$qbar) + 1)
  POOL[[i]]     <- c(pool$m, pool$qbar, pool$lower, pool$upper, pool$coverage)
}
AVG.cov.y0x <- Reduce("+", POOL) / length(POOL)
names(AVG.cov.y0x) <- c("m", "qbar", "lower", "upper", "coverage")
round(AVG.cov.y0x, 3)

# Evaluate cov of y1x
truth <- 0.5  
POOL <- list()
for (i in 1:nsim){
  #Extract means and variances
  Q            <- c(coefficient.var[[i]][[1]][2,3]/sqrt(coefficient.var[[i]][[1]][2,2] * coefficient.var[[i]][[1]][3,3]), 
                    coefficient.var[[i]][[2]][2,3]/sqrt(coefficient.var[[i]][[2]][2,2] * coefficient.var[[i]][[2]][3,3]),
                    coefficient.var[[i]][[3]][2,3]/sqrt(coefficient.var[[i]][[3]][2,2] * coefficient.var[[i]][[3]][3,3]), 
                    coefficient.var[[i]][[4]][2,3]/sqrt(coefficient.var[[i]][[4]][2,2] * coefficient.var[[i]][[4]][3,3]),
                    coefficient.var[[i]][[5]][2,3]/sqrt(coefficient.var[[i]][[5]][2,2] * coefficient.var[[i]][[5]][3,3])) %>% 
                    fisher.trans()
  
  U            <- rep(1 / (sample.size - 3), 5)
  #Pool the regular way
  pool         <- mice::pool.scalar(Q, U, n = 1000) # A really large number
  lower.bound  <- pool$qbar - qt(0.975, pool$df) * sqrt(pool$t)
  upper.bound  <- pool$qbar + qt(0.975, pool$df) * sqrt(pool$t)
  
  
  pool$lower   <- (exp(2 * lower.bound) - 1) / (exp(2 * lower.bound) + 1)
  pool$upper   <- (exp(2 * upper.bound) - 1) / (exp(2 * upper.bound) + 1)
  pool$coverage <- pool$lower < mean(truth) & mean(truth) < pool$upper
  pool$qbar <- (exp(2 * pool$qbar) - 1) / (exp(2 * pool$qbar) + 1)
  POOL[[i]]     <- c(pool$m, pool$qbar, pool$lower, pool$upper, pool$coverage)
}
AVG.cov.y1x <- Reduce("+", POOL) / length(POOL)
names(AVG.cov.y1x) <- c("m", "qbar", "lower", "upper", "coverage")
round(AVG.cov.y1x, 3)



OUT.mean <- (OUT[[1]]+OUT[[2]]+OUT[[3]]+OUT[[4]]+OUT[[5]])/5
OUT.mean <- Reduce('+', OUT)/5
hist(OUT.mean$y0[1:2500] - data$y0[1:2500], 
     xlab = "mean bias of missing y0", 
     main = "Histogram of mean bias of missing y0")
hist(OUT.mean$y1[2501:5000] - data$y1[2501:5000], 
     xlab = "mean bias of missing y1", 
     main = "Histogram of mean bias of missing y1")

c(OUT.mean$y0[1:2500] - data$y0[1:2500], OUT.mean$y1[2501:5000] - data$y1[2501:5000]) %>% hist(, xlab = "Mean bias", main = "Histogram of mean bias")

OUT.ind.test <- impute.multiple(20, incomplete.data, (0.8 -0.5*0.5) / (1 - 0.5^2))
OUT.ind.test.mean <- Reduce('+', OUT.ind.test)/20
c(OUT.ind.test.mean$y0[1:2500] - data$y0[1:2500], OUT.ind.test.mean$y1[2501:5000] - data$y1[2501:5000]) %>% var()
lm1<-lm(data$y0~data$y1+data$x)
lm1$residuals
var(lm1$residuals)