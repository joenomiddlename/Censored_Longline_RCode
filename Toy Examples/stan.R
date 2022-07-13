# stan.R - fit the censored method using Stan
source("simulate.R")
library(rstan)
library(brms)

# Need to create a variable with values 'none', left',
# 'right', or 'interval'

# We will use 'right' and 'none' here
y_BRMS <- c(y[1:100],
         ifelse(y[101:200] == 3, 2, y[101:200]),
         ifelse(y[201:300] == 2, 1, y[201:300]))
# NOTICE THAT FOR BRMS, WE NEED TO SUBTRACT 1 FROM THE
# RESPONSE AS THE CENSORSHIP INTERVAL IS CONSIDERED
# OPEN ON THE LEFT AND DATA HAVE TO BE INTEGERS

censored_txt <- rep('none', 300)
censored_txt[101:200] <- ifelse(y[101:200] == 3, 'right', 'none')
censored_txt[201:300] <- ifelse(y[201:300] == 2, 'right', 'none')

brms_mod <-
  brm(y_BRMS | cens(censored_txt) ~ x,
      data = data.frame(y_BRMS, x, censored_txt),
      family='Poisson')
summary(brms_mod)
# Compare the estimated intercept and linear value with the known values of a
#  and b.

# To make an interval censored value, simply define an upper
# bound. We use 10 here for demonstrative purposes.
censored_txt2 <- rep('none', 300)
censored_txt2[101:200] <- ifelse(y[101:200] == 3, 'interval', 'none')
censored_txt2[201:300] <- ifelse(y[201:300] == 2, 'interval', 'none')

# Define right interval endpoint as 10
U_BRMS2 <- y
U_BRMS2[censored_txt2 == 'interval'] <- 10

brms_mod2 <-
  brm(y_BRMS | cens(censored_txt2,U_BRMS2) ~ x,
      data = data.frame(y_BRMS, x, censored_txt2, U_BRMS2),
      family='Poisson')
summary(brms_mod2)
# Notice that the results are almost identical as before
