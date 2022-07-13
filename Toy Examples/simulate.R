# simulate.R - simulate some data
set.seed(29072021)
n <- 300
a <- 0
b <- 1
x <- rnorm(n, sd = 0.3)     # random observed covariates
lambda <- exp(a + b*x)
N <- rpois(n, lambda = lambda)
# Censor the data:
y <- c(N[1:100],
       ifelse(N[101:200]>=3, 3, N[101:200]),
       ifelse(N[201:300]>=2, 2, N[201:300]))
rm(N)  # To guarantee that the uncensored data not observed later
