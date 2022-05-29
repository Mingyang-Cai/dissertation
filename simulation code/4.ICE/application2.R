#application of ACTG175 prediction
rm(list=ls())
library(mice, warn.conflicts = FALSE)
library(miceadds)
library(MASS)
library(magrittr)
library(purrr)
library(speff2trial)



data(ACTG175)


testdata <- rbind(ACTG175[ACTG175$arms == 0,], ACTG175[ACTG175$arms == 3,])
testdata <- testdata[testdata$offtrt == 0, ]
testdata <- subset(testdata, select = -c(pidnum, zprior, strat, treat, 
                                         offtrt, days, cd496, r, cd820, cens))
testdata.0 <- testdata[testdata$arms == 0, ]
testdata.0$cd420.0 <- testdata.0$cd420
testdata.0$cd420.3 <- NA

testdata.3 <- testdata[testdata$arms == 3, ]
testdata.3$cd420.3 <- testdata.3$cd420
testdata.3$cd420.0 <- NA

testdata.oos <- rbind(ACTG175[ACTG175$arms == 1,], ACTG175[ACTG175$arms == 2,])
testdata.oos <- testdata.oos[testdata.oos$offtrt == 0, ]
testdata.oos <- subset(testdata.oos, select = -c(pidnum, zprior, strat, treat, 
                                                 offtrt, days, cd496, r, cd820, cens, arms, cd420))
testdata.oos$cd420.0 <- NA
testdata.oos$cd420.3 <- NA


incomplete.data <- rbind(testdata.0, testdata.3)
incomplete.data <- subset(incomplete.data, select = -c(arms, cd420))


covariate.balance <- matrix(NA, ncol = 4, nrow = ncol(incomplete.data) - 2)
for (i in 1 : nrow(covariate.balance)) {
  covariate.balance[i, 1] <- mean(incomplete.data[c(1: 316), i])
  covariate.balance[i, 2] <- sd(incomplete.data[c(1: 316), i])
  covariate.balance[i, 3] <- mean(incomplete.data[c(317: 693), i])
  covariate.balance[i, 4] <- sd(incomplete.data[c(317: 693), i])
}


#create a single imputed dataset
impute.single <- function(incomplete.data, prediction.data, rho.p){
  parm.y0 <- as.matrix(incomplete.data[, c(1:15)]) %>% 
    cbind(1, .) %>% 
    .norm.draw(incomplete.data$cd420.0, !is.na(incomplete.data$cd420.0), .)
  parm.y1 <- as.matrix(incomplete.data[, c(1:15)]) %>% 
    cbind(1, .) %>% 
    .norm.draw(incomplete.data$cd420.3, !is.na(incomplete.data$cd420.3), .)
  cov.y0y1.p <- rho.p * parm.y0$sigma * parm.y1$sigma
  
  
  mu.y0 <- as.matrix(prediction.data[, c(1:15)]) %>% 
    cbind(1, .) %>% multiply_by_matrix(parm.y0$beta)
  mu.y1 <- as.matrix(prediction.data[, c(1:15)]) %>% 
    cbind(1, .) %>% multiply_by_matrix(parm.y1$beta)
  
  
  
  for (i in 1 : nrow(prediction.data)) {
    temp <- mvrnorm(n = 1, mu = c(mu.y0[i], mu.y1[i]), Sigma = matrix(c(parm.y0$sigma^2, cov.y0y1.p, cov.y0y1.p, parm.y1$sigma^2), nrow = 2, ncol = 2))
    prediction.data$cd420.0[i] <- temp[1]
    prediction.data$cd420.3[i] <- temp[2]
    }
  
  
 
  return(prediction.data)
}
#create multiple imputed datasets
impute.multiple <- function(m = 5, incomplete.data, prediction.data, rho.p, ...){
  multi.imp.dat <- list()
  for (i in 1 : m) {
    multi.imp.dat[[i]] <-impute.single(incomplete.data, prediction.data, rho.p)  
  }
  return(multi.imp.dat)
}


OUT <- impute.multiple(5, incomplete.data, testdata.oos, 0.7)
OUT.mean <- Reduce('+', OUT)/5
ice.mean <- OUT.mean$cd420.3 - OUT.mean$cd420.0

plot(sort(ice.mean), col = "blue", pch = 0, ylim = c(-200, 250), cex = 0.6, xlab = "individual", ylab = "mean of ICE")



nrow(testdata.oos)


