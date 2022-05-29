#targeted learning
rm(list=ls())
library(mice, warn.conflicts = FALSE)
library(miceadds)
library(MASS)
library(magrittr)
library(purrr)
library(tmle)


rho.m <- 0.8                               #marginal correlation
#rho.p <- (0.8 -0.5*0.5) / (1 - 0.5^2)      #partial correlation
sample.size <- 5000
nsim <- 1000
set.seed(123)


y0.mean <- rep(0, nsim)
y1.mean <- rep(0, nsim)
y0.var <- rep(0, nsim)
y1.var <- rep(0, nsim)
y0y1.cor <- rep(0, nsim)
y0x.cor <- rep(0, nsim)
y1x.cor <- rep(0, nsim)

y0.mean.cr <- rep(0, nsim)
y1.mean.cr <- rep(0, nsim)
y0.var.cr <- rep(0, nsim)
y1.var.cr <- rep(0, nsim)
y0y1.cor.cr <- rep(0, nsim)
y0x.cor.cr <- rep(0, nsim)
y1x.cor.cr <- rep(0, nsim)




for (i in 1 : nsim) {
data <- mvrnorm(n = sample.size, mu = c(0, 1, 2), 
                Sigma = matrix(c(1, rho.m, 0.5, rho.m, 1, 0.5, 0.5, 0.5, 1), nrow = 3)) %>% as.data.frame
colnames(data)<- c("y0", "y1", "x")
incomplete.data <- data 
incomplete.data$y0[1 : 2500] <- NA
incomplete.data$y1[2501 : 5000] <- NA

#TMLE
Y <- c(incomplete.data$y1[1 : 2500], incomplete.data$y0[2501 : 5000])
A <- c(rep(1, 2500), rep(0, 2500))
W <- incomplete.data$x
result <- tmle(Y, A, W)
#summary(result)

y0.mean[i] <- mean(result$Qstar[,1])
y1.mean[i]<- mean(result$Qstar[,2])
y0.var[i] <- var(result$Qstar[,1])
y1.var[i] <- var(result$Qstar[,2])
y0y1.cor[i] <- cor(result$Qstar[,1], result$Qstar[,2])
y0x.cor[i] <- cor(result$Qstar[,1], W)
y1x.cor[i] <- cor(result$Qstar[,2], W)


y0.mean.boot <- rep(0, 100)
y1.mean.boot <- rep(0, 100)
y0.var.boot <- rep(0, 100)
y1.var.boot <- rep(0, 100)
y0y1.cor.boot <- rep(0, 100)
y0x.cor.boot <- rep(0, 100)
y1x.cor.boot <- rep(0, 100)
for (j in 1:100) {
  index <- sample(1:sample.size, replace = TRUE)
  data.boot <- data[index, ]
  incomplete.data <- data.boot 
  incomplete.data$y0[1 : 2500] <- NA
  incomplete.data$y1[2501 : 5000] <- NA
  
  Y <- c(incomplete.data$y1[1 : 2500], incomplete.data$y0[2501 : 5000])
  A <- c(rep(1, 2500), rep(0, 2500))
  W <- incomplete.data$x
  result <- tmle(Y, A, W)
  #summary(result)
  
  y0.mean.boot[j] <- mean(result$Qstar[,1])
  y1.mean.boot[j] <- mean(result$Qstar[,2])
  y0.var.boot[j] <- var(result$Qstar[,1])
  y1.var.boot[j] <- var(result$Qstar[,2])
  y0y1.cor.boot[j] <- cor(result$Qstar[,1], result$Qstar[,2])
  y0x.cor.boot[j] <- cor(result$Qstar[,1], W)
  y1x.cor.boot[j] <- cor(result$Qstar[,2], W)
}
y0.mean.boot <- y0.mean.boot[order(y0.mean.boot)]
y1.mean.boot <- y1.mean.boot[order(y1.mean.boot)]
y0.var.boot <- y0.var.boot[order(y0.var.boot)]
y1.var.boot <- y1.var.boot[order(y1.var.boot)]
y0y1.cor.boot <- y0y1.cor.boot[order(y0y1.cor.boot)]
y0x.cor.boot <- y0x.cor.boot[order(y0x.cor.boot)]
y1x.cor.boot <- y1x.cor.boot[order(y1x.cor.boot)]

if((y0.mean[i] > y0.mean.boot[2]) & (y0.mean[i] < y0.mean.boot[98])) y0.mean.cr[i] <- 1 
if((y1.mean[i] > y1.mean.boot[2]) & (y1.mean[i] < y1.mean.boot[98])) y1.mean.cr[i] <- 1 
if((y0.var[i] > y0.var.boot[2]) & (y0.var[i] < y0.var.boot[98])) y0.var.cr[i] <- 1 
if((y1.var[i] > y1.var.boot[2]) & (y1.var[i] < y1.var.boot[98])) y1.var.cr[i] <- 1
if((y0y1.cor[i] > y0y1.cor.boot[2]) & (y0y1.cor[i] < y0y1.cor.boot[98])) y0y1.cor.cr[i] <- 1
if((y0x.cor[i] > y0x.cor.boot[2]) & (y0x.cor[i] < y0x.cor.boot[98])) y0x.cor.cr[i] <- 1
if((y1x.cor[i] > y1x.cor.boot[2]) & (y1x.cor[i] < y1x.cor.boot[98])) y1x.cor.cr[i] <- 1
setTxtProgressBar(pb, i)
}
close(pb)


mean(y0.mean)
mean(y1.mean) 
mean(y0.var) 
mean(y1.var) 
mean(y0y1.cor) 
mean(y0x.cor) 
mean(y1x.cor)

sum(y0.mean.cr) / nsim
sum(y1.mean.cr) / nsim
sum(y0.var.cr) / nsim
sum(y1.var.cr) / nsim
sum(y0y1.cor.cr) / nsim 
sum(y0x.cor.cr) / nsim 
sum(y1x.cor.cr) / nsim 

