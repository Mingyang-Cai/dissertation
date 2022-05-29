#plot 
rm(list = ls())
library(mice, warn.conflicts = FALSE)
library(miceadds)
library(MASS)
library(magrittr)

sample.size <- 100
missingness <- 0.3

X <- rnorm(sample.size, mean = 0, sd = 1)
epsilon <- rnorm(sample.size, mean = 0, sd = 1)
Y <- X + X^2 + epsilon

complete.data <- data.frame(
  X,
  XX = X^2,
  Y
)

incomplete.data <- ampute(data=complete.data, prop=missingness, mech = "MAR", type = 'RIGHT', patterns = c(0, 0, 1),
                          weights = c(0, 0, 1))$amp

ini <- mice(incomplete.data, maxit = 0)
meth <- c("quadratic", "~I(X^2)", "")
pred <- ini$pred
pred[, "XX"] <- 0

# Impute data
imp.pc <- mice(incomplete.data, meth = meth, pred = pred, quad.outcome = "Y", print = FALSE)

#complete.data$X, complete.data$Y  gray50
#incomplete.data[!is.na(incomplete.data$X), ]$X, incomplete.data[!is.na(incomplete.data$X), ]$Y col="#006CC2B3"
#complete(imp.pc)$X, complete(imp.pc)$Y  B61A51B3

plot(complete.data$X, complete.data$Y, xlab = "X", ylab = "Y", col = "gray50", pch = 16)
points(incomplete.data[!is.na(incomplete.data$X), ]$X, incomplete.data[!is.na(incomplete.data$X), ]$Y, 
       col = "#006CC2B3", pch = 16)
legend(x="topleft", legend = c("complete", "observed"), col = c("gray50", "#006CC2B3"), pch = c(16, 16))




plot(complete.data$X, complete.data$Y, xlab = "X", ylab = "Y", col = "gray50", pch = 16)
points(complete(imp.pc)$X, complete(imp.pc)$Y, col = "#B61A51B3", pch = 16)
legend(x="topleft", legend = c("complete", "imputed"), col = c("gray50", "#B61A51B3"), pch = c(16, 16))



plot(density(incomplete.data[!is.na(incomplete.data$X), ]$X), col = "#006CC2B3", xlab = "X", main = "Density of X")
lines(density(complete.data$X), col = "gray50")
legend(x="topleft", legend = c("complete", "observed"), col = c("gray50", "#006CC2B3"), lty = c(1, 1))

plot(density(complete(imp.pc)$X), col = "#B61A51B3", xlab = "X", main = "Density of X")
lines(density(complete.data$X), col = "gray50")
legend(x="topleft", legend = c("complete", "imputed"), col = c("gray50", "#B61A51B3"), lty = c(1, 1))
