# r-inla.R - fit the censored method using R-INLA
source("simulate.R")
library(INLA)
#  version 22.04.16 or higher.

# Define the lower bounds for the 300 censorship intervals.
# Note that for values of low set to Inf, INLA assumes the
# data are perfectly observed (no censorship).
low <- c(rep(Inf, 100),
         ifelse(y[101:200] == 3, 3, Inf),
         ifelse(y[201:300] == 2, 2, Inf))

# Define the upper bounds, which are always infinite.
high <- rep(Inf, n)

# What do this data look like?
cbind(low, y, high)

# Fit the model
INLA_mod <- inla(inla.mdata(cbind(y, low, high)) ~ 1 + x,
                 family = "cenpoisson2",
                 data = data.frame(y, low, high, x))

# View the results
summary(INLA_mod)

# If interval-censoring is desired change the upper bound
# Suppose the upper bound is 10
high <- c(y[1:100],
          ifelse(y[101:200] == 3, 10, y[101:200]),
          ifelse(y[201:300] == 2, 10, y[201:300]))

# Refit the model
INLA_mod <- inla(inla.mdata(cbind(y, low, high)) ~ 1 + x,
                 family = "cenpoisson2",
                 data = data.frame(y, low, high, x))

# View the results
summary(INLA_mod)
# Doesn't change much
