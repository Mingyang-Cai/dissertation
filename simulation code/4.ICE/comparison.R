#comparision new
rm(list=ls())
library(mice, warn.conflicts = FALSE)
library(miceadds)
library(MASS)
library(magrittr)
library(purrr)
library(tmle)

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



rho.m <- 0.8                               #marginal correlation
#rho.p <- (0.8 -0.5*0.5) / (1 - 0.5^2)      #partial correlation
sample.size <- 5000
index.y0 <- c(1 : 25) * 100
index.y1 <- c(26 : 50) * 100
true.var <- (1-0.73^2)
set.seed(123)

data <- mvrnorm(n = sample.size, mu = c(0, 1, 2), 
                Sigma = matrix(c(1, rho.m, 0.5, rho.m, 1, 0.5, 0.5, 0.5, 1), nrow = 3)) %>% as.data.frame
colnames(data)<- c("y0", "y1", "x")
true.ice <- data$y1 - data$y0
incomplete.data <- data 
incomplete.data$y0[1 : 2500] <- NA
incomplete.data$y1[2501 : 5000] <- NA
#FCS approach rho.p = 0
OUT.1 <- impute.multiple(20, incomplete.data, 0)
OUT.1.mean <- Reduce('+', OUT.1)/20
est.ice.1 <- OUT.1.mean$y1 - OUT.1.mean$y0
bias.1 <- est.ice.1 - true.ice
hist(bias.1, xlab = "mean of bias", main = " ")
mean(bias.1)
var(bias.1)


plot(data$y1[index.y0], data$y1[index.y0] - data$y0[index.y0], col = "#006CC2B3", 
     pch = 18, xlim = c(-2, 3) , ylim = c(-3, 5), xlab = (expression(Y[obs])), ylab = "ICE")
points(data$y0[index.y1], data$y1[index.y1] - data$y0[index.y1], col = "#006CC2B3", pch = 20)
for (i in 1 : 20) {
  points(data$y1[index.y0], data$y1[index.y0] - OUT.1[[i]]$y0[index.y0], col = "#B61A51B3", pch = 18)
  points(data$y0[index.y1], OUT.1[[i]]$y1[index.y1] - data$y0[index.y1], col = "#B61A51B3", pch = 20)
}
legend("topleft", c("control", "treatment"), pch = c(20, 18))


OUT.1.index.data <- matrix(NA, 20, 50)
for (i in 1 : 20) {
  OUT.1.index.data[i, ] <- c(data$y1[index.y0] - OUT.1[[i]]$y0[index.y0], OUT.1[[i]]$y1[index.y1] - data$y0[index.y1])
}
apply(OUT.1.index.data, 2, var) 


#FCS approach rho.p = 0.73
OUT.2 <- impute.multiple(20, incomplete.data, 0.73)
OUT.2.mean <- Reduce('+', OUT.2)/20
est.ice.2 <- OUT.2.mean$y1 - OUT.2.mean$y0
bias.2 <- est.ice.2 - true.ice
hist(bias.2, xlab = "mean of bias", main = " ")
mean(bias.2)
var(bias.2)
plot(data$y1[index.y0], data$y1[index.y0] - data$y0[index.y0], col = "#006CC2B3", 
     pch = 18, xlim = c(-2, 3) , ylim = c(-3, 5), xlab = (expression(Y[obs])), ylab = "ICE")
points(data$y0[index.y1], data$y1[index.y1] - data$y0[index.y1], col = "#006CC2B3", pch = 20)
for (i in 1 : 20) {
  points(data$y1[index.y0], data$y1[index.y0] - OUT.2[[i]]$y0[index.y0], col = "#B61A51B3", pch = 18)
  points(data$y0[index.y1], OUT.2[[i]]$y1[index.y1] - data$y0[index.y1], col = "#B61A51B3", pch = 20)
}
legend("topleft", c("control", "treatment"), pch = c(20, 18))
OUT.2.index.data <- matrix(NA, 20, 50)
for (i in 1 : 20) {
  OUT.2.index.data[i, ] <- c(data$y1[index.y0] - OUT.2[[i]]$y0[index.y0], OUT.2[[i]]$y1[index.y1] - data$y0[index.y1])
}
apply(OUT.2.index.data, 2, var)%>%mean


#FCS approach rho.p = 0.99
OUT.3 <- impute.multiple(20, incomplete.data, 0.99)
OUT.3.mean <- Reduce('+', OUT.3)/20
est.ice.3 <- OUT.3.mean$y1 - OUT.3.mean$y0
bias.3 <- est.ice.3 - true.ice
hist(bias.3, xlab = "mean of bias", main = " ")
mean(bias.3)
var(bias.3)
plot(data$y1[index.y0], data$y1[index.y0] - data$y0[index.y0], col = "#006CC2B3", 
     pch = 18, xlim = c(-2, 3) , ylim = c(-3, 5), xlab = (expression(Y[obs])), ylab = "ICE")
points(data$y0[index.y1], data$y1[index.y1] - data$y0[index.y1], col = "#006CC2B3", pch = 20)
for (i in 1 : 20) {
  points(data$y1[index.y0], data$y1[index.y0] - OUT.3[[i]]$y0[index.y0], col = "#B61A51B3", pch = 18)
  points(data$y0[index.y1], OUT.3[[i]]$y1[index.y1] - data$y0[index.y1], col = "#B61A51B3", pch = 20)
}
legend("topleft", c("control", "treatment"), pch = c(20, 18))
OUT.3.index.data <- matrix(NA, 20, 50)
for (i in 1 : 20) {
  OUT.3.index.data[i, ] <- c(data$y1[index.y0] - OUT.3[[i]]$y0[index.y0], OUT.3[[i]]$y1[index.y1] - data$y0[index.y1])
}
apply(OUT.3.index.data, 2, var)%>%mean

#TMLE
Y <- c(incomplete.data$y1[1 : 2500], incomplete.data$y0[2501 : 5000])
A <- c(rep(1, 2500), rep(0, 2500))
W <- as.matrix(incomplete.data$x)
OUT.4 <- list()
for (i in 1 : 20) {
  result <- tmle(Y, A, W)
  OUT.4[[i]] <- cbind(result$Qstar[,2], result$Qstar[,1])
}
OUT.4.mean <- Reduce('+', OUT.4)/20
est.ice.4 <- OUT.4.mean[, 1] - OUT.4.mean[, 2]
bias.4 <- est.ice.4 - true.ice
hist(bias.4, xlab = "mean of bias", main = " ")
mean(bias.4)
var(bias.4)
plot(data$y1[index.y0], data$y1[index.y0] - data$y0[index.y0], col = "#006CC2B3", 
     pch = 18, xlim = c(-2, 3) , ylim = c(-3, 5), xlab = (expression(Y[obs])), ylab = "ICE")
points(data$y0[index.y1], data$y1[index.y1] - data$y0[index.y1], col = "#006CC2B3", pch = 20)
for (i in 1 : 20) {
  points(data$y1[index.y0], data$y1[index.y0] - OUT.4[[i]][, 1][index.y0], col = "#B61A51B3", pch = 18)
  points(data$y0[index.y1], OUT.4[[i]][,1][index.y1] - data$y0[index.y1], col = "#B61A51B3", pch = 20)
}
legend("topleft", c("control", "treatment"), pch = c(20, 18))
OUT.4.index.data <- matrix(NA, 20, 50)
for (i in 1 : 20) {
  OUT.4.index.data[i, ] <- c(data$y1[index.y0] - OUT.4[[i]][,2][index.y0], OUT.4[[i]][,1][index.y1] - data$y0[index.y1])
}
apply(OUT.4.index.data, 2, var)%>%mean
