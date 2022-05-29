#application of ACTG175  imputation
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
                                        offtrt, cd420, cd496, r, cd820, cens))
testdata.0 <- testdata[testdata$arms == 0, ]
testdata.0$days.0 <- testdata.0$days
testdata.0$days.3 <- NA

testdata.3 <- testdata[testdata$arms == 3, ]
testdata.3$days.3 <- testdata.3$days
testdata.3$days.0 <- NA

incomplete.data <- rbind(testdata.0, testdata.3)
incomplete.data <- subset(incomplete.data, select = -c(arms, days))


covariate.balance <- matrix(NA, ncol = 4, nrow = ncol(incomplete.data) - 2)
for (i in 1 : nrow(covariate.balance)) {
  covariate.balance[i, 1] <- mean(incomplete.data[c(1: 316), i])
  covariate.balance[i, 2] <- sd(incomplete.data[c(1: 316), i])
  covariate.balance[i, 3] <- mean(incomplete.data[c(317: 693), i])
  covariate.balance[i, 4] <- sd(incomplete.data[c(317: 693), i])
}


#create a single imputed dataset
impute.single <- function(incomplete.data, rho.p){
  parm.y0 <- as.matrix(incomplete.data[, c(1:15)]) %>% 
    cbind(1, .) %>% 
    .norm.draw(incomplete.data$days.0, !is.na(incomplete.data$days.0), .)
  parm.y1 <- as.matrix(incomplete.data[, c(1:15)]) %>% 
    cbind(1, .) %>% 
    .norm.draw(incomplete.data$days.3, !is.na(incomplete.data$days.3), .)
  cov.y0y1.p <- rho.p * parm.y0$sigma * parm.y1$sigma
  #impute y0
  x0 <- as.matrix(incomplete.data[, c(1:15)]) %>% 
    cbind(1, .) %>% 
    extract(is.na(incomplete.data$days.0), )
  imp.y0 <- x0 %>% 
    multiply_by_matrix(parm.y0$beta) %>%
    add((as.matrix(incomplete.data$days.3[is.na(incomplete.data$days.0)]) - x0 %*% parm.y1$beta) * cov.y0y1.p / (parm.y1$sigma)^2) %>%
    add(rnorm(sum(is.na(incomplete.data$days.0))) * sqrt((1 - rho.p^2)) * parm.y0$sigma)
  #impute y1
  x1 <- as.matrix(incomplete.data[, c(1:15)]) %>% 
    cbind(1, .) %>% 
    extract(is.na(incomplete.data$days.3), )
  imp.y1 <- x1 %>% 
    multiply_by_matrix(parm.y1$beta) %>%
    add((as.matrix(incomplete.data$days.0[is.na(incomplete.data$days.3)]) - x1 %*% parm.y0$beta) * cov.y0y1.p / (parm.y0$sigma)^2) %>%
    add(rnorm(sum(is.na(incomplete.data$days.3))) * sqrt((1 - rho.p^2)) * parm.y1$sigma)
  incomplete.data$days.0[317 : 693] <- imp.y0
  incomplete.data$days.3[1 : 316] <- imp.y1
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


OUT <- impute.multiple(5, incomplete.data, 0.7)
OUT.mean <- Reduce('+', OUT)/5
ice.mean <- OUT.mean$days.3 - OUT.mean$days.0
ice.mean.cg <- sort(ice.mean[1 : 316])
ice.mean.tg <- sort(ice.mean[317 : 693])
plot(ice.mean.cg, col = "blue", pch = 0, cex = 0.6, ylim = c(-1000, 1000), xlab = "individual", ylab = "mean of ICE")
points(ice.mean.tg, col = "red", pch = 6, cex = 0.6)
legend("topleft", c("treatment A", "treatment B"), col = c("blue", "red"), pch = c(0, 6))
