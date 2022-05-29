#add plot with stef's comment
rm(list=ls())
library(mice, warn.conflicts = FALSE)
library(miceadds)
library(MASS)
library(magrittr)
library(purrr)
library(tmle)
library(ggplot2)

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


#FCS approach rho.p = 0.73
OUT.2 <- impute.multiple(20, incomplete.data, 0.73)
OUT.2.mean <- Reduce('+', OUT.2)/20
est.ice.2 <- OUT.2.mean$y1 - OUT.2.mean$y0
est.ice.2.cg <- est.ice.2[2501 : 5000]
est.ice.2.tg <- est.ice.2[1 : 2500]
true.ice <- data$y1 - data$y0


data.plot <- data.frame(bias = est.ice.2,
                        obs = c(data$y1[1:2500], data$y0[2501:5000]),
                        group = c(rep("treatment", 2500), rep("control", 2500)))

ggplot(data.plot, aes(obs, bias, colour = group)) + geom_smooth(method = "loess") + xlab(expression(Y[obs]))+ylab("Y(1) - Y(0)")





cg <- data.frame(
  bias = est.ice.2.cg - true.ice[2501 : 5000],
  obs = data$y0[2501 : 5000]
)

tg <- data.frame(
  bias = est.ice.2.tg - true.ice[1 : 2500],
  obs = data$y1[1 : 2500]
)
  
  
plot(cg$obs, cg$ice)
cg[order(cg$obs),]$ice
ggplot() + geom_smooth(data = cg, aes(obs, bias), method = "loess", color = "blue") +geom_smooth(data = tg, aes(obs, bias), method = "loess",  color = "red") + xlab(expression(Y[obs]))+
  scale_color_discrete(name = "the colour")








bias.2 <- est.ice.2 - true.ice
hist(bias.2, xlab = "mean of bias", main = " ")
mean(bias.2)
var(bias.2)
plot(c(data$y1[index.y0] - OUT.2[[1]]$y0[index.y0], OUT.2[[1]]$y1[index.y1] - data$y0[index.y1]) , ylim = c(-3, 5), col = "#B61A51B3", pch = 20, xlab = " ", ylab = " ")
for (i in 2 : 20) {
  points(c(data$y1[index.y0] - OUT.2[[i]]$y0[index.y0], OUT.2[[i]]$y1[index.y1] - data$y0[index.y1]), col = "#B61A51B3", pch = 20)
}
points(c(data$y1[index.y0] - data$y0[index.y0], data$y1[index.y1] - data$y0[index.y1]), col = "#006CC2B3")
OUT.2.index.data <- matrix(NA, 20, 50)
for (i in 1 : 20) {
  OUT.2.index.data[i, ] <- c(data$y1[index.y0] - OUT.2[[i]]$y0[index.y0], OUT.2[[i]]$y1[index.y1] - data$y0[index.y1])
}
apply(OUT.2.index.data, 2, var)%>%mean

aa <- data$y0 - OUT.2[[1]]$y0
aa[order(data$y0),]



plot(NULL, xlim = c(0, 5000), ylim = c(-1, 1), ylab = "bias", xlab = "index")
lines(lowess(sort(data$y0 - OUT.2[[1]]$y0)), col = "red")
lines(lowess(sort(OUT.2[[1]]$y1 - data$y1, decreasing = TRUE)), col = "blue")
legend(1500, 0.75, c("control", "treatment"), col = c("blue","red"), lty = c(1, 1))

aa <- data.frame(
  ice = data$y0 - OUT.2[[1]]$y0,
  obs = data$y0
)



bb <- OUT.2[[1]]$y1 - data$y1
#bb[order(data$y1, decreasing = TRUE)]
plot(NULL, xlim = c(-5, 10), ylim = c(-1, 1), ylab = "bias", xlab = "index")
lines(lowess(aa[order(aa$obs), ]$ice), col = "red")
lines(lowess(bb[order(data$y1, decreasing = TRUE)]), col = "blue")
legend(1500, 0.75, c("control", "treatment"), col = c("blue","red"), lty = c(1, 1))

