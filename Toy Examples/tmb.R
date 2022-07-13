# tmb.R - fit the censored method using TMB
source("simulate.R")
library(sdmTMB)
# Install from https://github.com/pbs-assess/sdmTMB

# Define the lower and upper bounds. Note that equality
# in bounds is equivalent to no censoring
low <- y

high <- c(y[1:100],
         ifelse(y[101:200] == 3, NA, y[101:200]),
         ifelse(y[201:300] == 2, NA, y[201:300]))

# Fit the model (note that an arbitrary mesh needs to be made)
# The mesh is not used in the computation with spatial='off'
sdmTMB_mod <- sdmTMB(
  data = data.frame(x=x,y=y),
  formula = y ~ x,
  mesh = make_mesh(data.frame(x=1:300,y=1:300),
                   xy_cols = c('x','y'),
                   n_knots = 3),
  family = censored_poisson(link = "log"),
  experimental = list(lwr = low, upr = high),
  spatial = 'off',
  spatiotemporal = 'off'
)

# View the results
summary(sdmTMB_mod)

# If interval-censoring is desired change the upper bound
high <- c(y[1:100],
          ifelse(y[101:200] == 3, 10, y[101:200]),
          ifelse(y[201:300] == 2, 10, y[201:300]))

# Refit the model
sdmTMB_mod <- sdmTMB(
  data = data.frame(x=x,y=y),
  formula = y ~ x,
  mesh = make_mesh(data.frame(x=1:300,y=1:300),
                   xy_cols = c('x','y'),
                   n_knots = 3),
  family = censored_poisson(link = "log"),
  experimental = list(lwr = low, upr = high),
  spatial = 'off',
  spatiotemporal = 'off'
)

# View the results
summary(sdmTMB_mod)
# Doesn't change much
