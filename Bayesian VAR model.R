library(bvartools)
library(readr)
library(fpp2)
library(tseries)
library(ggplot2)
library(forecast)
library(zoo)
library(xts)
library(vars)
library(urca)
library(sandwich)
# Dropping the date column so I can declare data as time series
dg = subset(Complete_Quarterly_Dataset, select=-c(Date))
ds <- ts(dg, start = c(1960,2), end=c(2021,2), frequency = 4)
# This function produces an object, which contains information on the specification of the VAR model that should be estimated.
model <- gen_var(ds, p=1, deterministic = "const",
                 iterations = 5000, burnin = 1000)
# Adding Model Priors
model_with_priors <- add_priors(model,
                                coef = list(v_i = 0, v_i_det = 0),
                                sigma = list(df = 1, scale = .0001))
# The output of add_priors can be used as the input for user-written algorithms for posterior simulation.
# However, bvar comes with built-in posterior simulation functions.
bvar_est <- draw_posterior(model_with_priors)
# Setting up a Gibbs sampler algorithm 
# Reset random number generator for reproducibility 
set.seed(1234567)

y <- t(model_with_priors$data$Y) # Matrix of dependent variables
x <- t(model_with_priors$data$Z) # Matrix of regressors for the model y_t = Ax_t + u_t

tt <- ncol(y) # Number of observations
k <- nrow(y) # Number of endogenous variables
m <- k * nrow(x) # Number of estimated coefficients

# Set (uninformative) priors
a_mu_prior <- model_with_priors$priors$coefficients$mu # Vector of prior parameter means
a_v_i_prior <- model_with_priors$priors$coefficients$v_i # Inverse of the prior covariance matrix

u_sigma_df_prior <- model_with_priors$priors$sigma$df # Prior degrees of freedom
u_sigma_scale_prior <- model_with_priors$priors$sigma$scale # Prior covariance matrix
u_sigma_df_post <- tt + u_sigma_df_prior # Posterior degrees of freedom

# Initial values
u_sigma <- diag(.00001, k)
u_sigma_i <- solve(u_sigma)

# Data containers for posterior draws
draws_a <- matrix(NA, m, iterations)
draws_sigma <- matrix(NA, k^2, iterations)

# Number of iterations of the Gibbs sampler
iterations <- model_with_priors$model$iterations
# Number of burn-in draws
burnin <- model_with_priors$model$burnin
# Total number of draws
draws <- iterations + burnin

# Storage for posterior draws
draws_a <- matrix(NA, m, iterations)
draws_sigma <- matrix(NA, k^2, iterations)

# Start Gibbs sampler
for (draw in 1:draws) {
  # Draw conditional mean parameters
  a <- post_normal(y, x, u_sigma_i, a_mu_prior, a_v_i_prior)
  
  # Draw variance-covariance matrix
  u <- y - matrix(a, k) %*% x # Obtaining residuals
  u_sigma_scale_post <- solve(u_sigma_scale_prior + tcrossprod(u))
  u_sigma_i <- matrix(rWishart(1, u_sigma_df_post, u_sigma_scale_post)[,,1], k)
  u_sigma <- solve(u_sigma_i) # Invert Sigma_i to obtain sigma
  
  # Store draws
  if (draw > burnin) {
    draws_a[, draw - burnin] <- a
    draws_sigma[, draw - burnin] <- solve(u_sigma_i) # Invert sigma_i to obtain sigma
  }
}
# Collecting relevant output of Gibbs sampler in a standardized object
bvar_est_two <- bvar(y = model_with_priors$data$Y,
                 x = model_with_priors$data$Z,
                 A = draws_a[1:100,],
                 C = draws_a[101:110,],
                 Sigma = draws_sigma)

# Summary statistics
summary(bvar_est)
# Posterior draws can be thinned with the function thin:
bvar_est <- thin(bvar_est, thin = 10)
summary(bvar_est)
# Making predictions
bvar_pred <- predict(bvar_est, n.ahead = 4)
plot(bvar_pred)

# Experiencing difficulty with the code below.
plot(bvar_est)
plot(bvar_est, type = "trace")
