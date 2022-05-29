#section 3 mpmm  x1+x2>=3
rm(list=ls())
library(mice, warn.conflicts = FALSE)
library(miceadds)
library(MASS)
library(magrittr)

set.seed(123)
sample.size <- 2000
missingness <- 0.3
simulations <- 1000
#correlation.coef <- 0.8

#function for pmm-cra
mice.impute.pmm.cra <- function(data, format = "imputes", ...){
  order <- dimnames(data)[[1]]
  res <- pmm.cra.impute(data, ...)
  return(mice:::single2imputes(res[order,], is.na(data)))
}
pmm.cra.impute<-function(data, ...){
  data <- as.data.frame(data)
  r <- 1 * is.na(data)
  pat <- apply(r, 1, function(x) paste(as.numeric(x), collapse=''))
  r <- matrix(r[order(pat), ], dim(data))
  data <- data[order(pat), , drop = FALSE]
  pat <- pat[order(pat)]
  nmpat <- length(unique(pat)) - 1
  pat.index <- c(match(unique(pat), pat), nrow(data) + 1)[-1]
  #impute for different missing patterns
  for (i in 1 : nmpat) {                                  
    dataa <- data[c(1 : (pat.index[i + 1] - 1)), , drop = FALSE]
    rr <- unique(r)[i + 1, ]
    xx <- dataa[, which(rr == 0), drop = FALSE]
    yy <- dataa[, which(rr == 1), drop = FALSE]
    ry <- !is.na(yy)[, 1]
    data[c(pat.index[i] : (pat.index[i + 1] - 1)), ] <- pmm.cra.single(yy, ry, xx, dimnames(data)[[2]])
  }
  return(data)
}
pmm.cra.single <- function(y, ry, x, col.order, wy = NULL, ...){
  if (is.null(wy)) wy <- !ry
  ES <- eigen(solve(cov(y[ry, ], y[ry, ])) %*% cov(y[ry, ], x[ry,]) 
              %*% solve(cov(x[ry,], x[ry,])) %*% cov(x[ry,], y[ry, ]))
  alpha <- ES$vectors[,1]
  alpha <- as.matrix(alpha)
  frv.obs <- as.matrix(y[ry, ]) %*% alpha
  frv.mis <- as.matrix(y[wy, ]) %*% alpha
  frv.total <- rbind(frv.obs, frv.mis)
  imp<-mice.impute.pmm(frv.total,!is.na(frv.total),x)
  frv.mis.est <- as.matrix(imp)
  y[wy, ] <- y[ry, ,drop = FALSE][match(frv.mis.est, frv.obs),]
  output <- cbind(y, x)
  return(output[wy, col.order])
}

#data generation function
create.data <- function(n) {
  mu.x <- c(0, 1)
  sigma.x <- matrix(data = c(4, 3.2, 3.2, 4), nrow = 2, byrow = T)
  X <- mvrnorm(n, mu = mu.x, Sigma = sigma.x)
  X1 <- X[, 1]
  X3 <- X[, 2]
  X2 <- 3 - X1 + runif(n, 0, 1)
  return(cbind(X1, X2, X3))
}

#mean(X1) = 0
#mean(X2) = 3.5

#output project
OUT <- list()
#simulation
pb <- txtProgressBar(min = 0, max = simulations, style = 3)
for (i in 1 : simulations) {
  #data generation
  complete.data <- create.data(sample.size)
  #ampute
  #incomplete.data <- ampute(complete.data, prop = missingness, mech = "MCAR", 
  #                         patterns = c(0, 0, 1))$amp
  incomplete.data <- ampute(complete.data, prop = missingness, mech = "MAR", 
                            type = 'RIGHT', patterns = c(0, 0, 1),
                            weights = c(0, 0, 1))$amp
  colnames(incomplete.data) <- c("X1", "X2", "X3")
  #impute
  blk <- list(c("X1", "X2"), c("X3"))
  method<- c("pmm.cra", "")
  OUT[[i]] <- mice(incomplete.data, m = 5, blocks = blk, method = method, print = FALSE)
  setTxtProgressBar(pb, i)
}
close(pb)

# Evaluate E(X1)
evaluate.sims <- function(sims, truth = 0){
  POOL <- list()
  pb <- txtProgressBar(min = 0, max = simulations, style = 3)
  for (i in 1:length(sims)){
    #Extract means and variances
    Q            <- unlist(with(sims[[i]], mean(X1))$analyses)
    U            <- unlist(with(sims[[i]], var(X1))$analyses) / nrow(sims[[i]]$data)
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
EVAL <- evaluate.sims(OUT)

# Summarize
AVG.EVAL <- Reduce("+", EVAL) / length(EVAL)
round(AVG.EVAL, 3)

# Evaluate E(X2)
evaluate.sims <- function(sims, truth = 3.5){
  POOL <- list()
  pb <- txtProgressBar(min = 0, max = simulations, style = 3)
  for (i in 1:length(sims)){
    #Extract means and variances
    Q            <- unlist(with(sims[[i]], mean(X2))$analyses)
    U            <- unlist(with(sims[[i]], var(X2))$analyses) / nrow(sims[[i]]$data)
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
EVAL <- evaluate.sims(OUT)

# Summarize
AVG.EVAL <- Reduce("+", EVAL) / length(EVAL)
round(AVG.EVAL, 3)


test <-OUT %>% map(~.x %>% complete("all") %>% 
                     map(~.x %>% extract(c("X1", "X2"))%>% rowSums() %>% 
                           is_greater_than(3)%>%sum())) %>% unlist() %>% mean()



test