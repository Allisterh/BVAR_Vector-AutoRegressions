install.packages("shrinkTVP")
install.packages("RhpcBLASctl")
library(RhpcBLASctl)
library(foreach)
# Running the model
data("usmacro.update")
library(shrinkTVP)
set.seed(123)
sim <- simTVP(theta = c(0.2,0,0), beta_mean = c(-1.5, -0.3, 0))
data <- sim$data
res <- shrinkTVP(y ~ x1 +x2, data = Complete_Quarterly_Dataset)
# Specifying Priors
# Fixing the pole parameters
res_hierlasso <- shrinkTVP(y ~ x1 + x2, data = data,
                           learn_a_xi = FALSE, learn_a_tau = FALSE,
                           a_xi = 1, a_tau = 1, display_progress = FALSE)
# Fixing the global shrinkage parameters
res_lasso <- shrinkTVP(y ~ x1 + x2, data = data,
                       learn_a_xi = FALSE, learn_a_tau = FALSE, a_xi = 1, a_tau = 1,
                       learn_kappa2_B = FALSE, learn_lambda2_B = FALSE,
                       kappa2_B = 100, lambda2_B = 100,
                       display_progress = FALSE)
# Change the prior type
res_tg <- shrinkTVP(y ~ x1+ x2, data = data,
                    mod_type = "triple",
                    display_progress = FALSE)
# Fixing the tail parameters
res_hs <- shrinkTVP(y ~ x1 + x2, data = data,
                    mod_type = "triple",
                    learn_a_xi = FALSE, learn_a_tau = FALSE, a_xi = 0.5, a_tau = 0.5,
                    learn_c_xi = FALSE, learn_c_tau = FALSE, c_xi = 0.5, c_tau = 0.5,
                    learn_kappa2_B = FALSE, learn_lambda2_B = FALSE,
                    display_progress = FALSE)
# Stochastic volatility specification
res_sv <- shrinkTVP(y ~ x1 + x2, data = data, sv = TRUE,
                    display_progress = FALSE)
# Specifying the hyperparameters
res_hyp <- shrinkTVP(y ~ x1 + x2, data = data,
                     hyperprior_param = list(beta_a_xi = 5, alpha_a_xi = 10),
                     display_progress = FALSE)
# Tuning the Metropolis-Hastings steps
res_MH <- shrinkTVP(y ~ x1 + x2, data = data,
                    MH_tuning = list(a_xi_adaptive = FALSE,
                                     a_tau_max_adapt = 0.001,
                                     a_tau_batch_size = 20),
                    display_progress = FALSE)
# Posterior inference: Summarize and vizualise the posterior distribtuion
summary(res, showprior = FALSE)
plot(res)
# Predictive performances and model comparison
library("RColorBrewer")
color <- brewer.pal(5, "RdBu")
plot(res, pars = "beta", xlim = c(100,200),
     probs =  seq(0.1, 0.9, by = 0.1),
     quantiles = T, quantcol = color[5:2], quantlty = 1,
     quantlwd = 3, col = color[1], lwd = 3, shadecol = "gold1")
plot(res, pars = "theta_sr")
res_LPDS <- shrinkTVP(y ~ x1 + x2, data = data[1:(nrow(data) -1),],
                      display_progress = FALSE)
LPDS(res_LPDS, data[nrow(data),])
eval_pred_dens(1:3, res_LPDS, data[nrow(data),])
curve(eval_pred_dens(x, res_LPDS, data[nrow(data),]), to = 12,
      ylab = bquote("p("* y[t[0]+1] * "\uff5c" * y[1] * ","
                    * ldots * "," ~ y[t[0]] * "," ~ x[t[0]+1] * ")"),
      xlab = expression(y[t[0]+1]), lwd = 2.5, col = "skyblue", axes = FALSE)
abline(v = data$y[nrow(data)])
axis(1)
axis(2)
# Predictive exercise on my Dataset
install.packages("bvarsv")
library("bvarsv")
dsg = subset(Complete_Quarterly_Dataset, select = -c(Date))
mym <- ts(us_data, start = c(1960,2), end=c(2021,2), frequency = 4)

lags <- dsg[1:(nrow(Complete_Quarterly_Dataset)-1),]
colnames(lags) <- paste0(colnames(lags), "_lag")
us_data <- data.frame(y1 = dsg[2:nrow(dsg), "y1"],
                      lags)
us_res <- shrinkTVP(y1 ~ y1_lag + y2_lag + y3_lag + y4_lag + y5_lag + y6_lag + y7_lag + y8_lag + y9_lag + y10_lag, us_data,
                    niter = 60000, nburn = 10000, nthin = 10,
                    display_progress = FALSE)
summary(us_res, showprior = FALSE)

# Calculating LPDS in Multicore
library(doParallel)
library(foreach)
library(zoo)
library(RhpcBLASctl)
install.packages("doSNOW")
install.packages("doParallel")
install.packages("Rmpi")
install.packages("doMPI")
library(doSNOW)
library(doMPI)
library(foreach)
# Defining how many periods to calculate LPDS for
Tmax <- nrow(us_data) - 1
T0 <- Tmax - 44

# Determine number of cores to be used and register parallel backend 
ncores <- 4
cl <- makeCluster(ncores)
registerDoParallel(cl)

lpds <- foreach(t = T0:Tmax, .combine = "cbind",
                .packages = c("RhpcBLASctl", "shrinkTVP"),
                .errorhandling = "pass") %do% {
                  
                  set.seed(t)
                  
                  niter <- 30000
                  nburn <- 15000
                  nthin <- 5
                  
                  # Set number of BLAS threads, so they dont interfere with each other
                  blas_set_num_threads(1)
                  
                  # Create data_t from all data up to time t and 
                  # y_test and x_test from data at time t+1
                  data_test <- us_data[t+1,]
                  data_t <- us_data[1:t,]
                  
                  #Run MCMC to calculate all LPDS
                  # Fully hierarchical triple gamma
                  res_FH_TG <- shrinkTVP(y1 ~ y1_lag + y2_lag + y3_lag + y4_lag + y5_lag + y6_lag + y7_lag + y8_lag + y9_lag + y10_lag, data = data_t,
                                         niter = niter, nburn = nburn, nthin = nthin)
                  
                  # Hierarchical triple gamma
                  res_H_TG <- shrinkTVP(y1 ~ y1_lag + y2_lag + y3_lag + y4_lag + y5_lag + y6_lag + y7_lag + y8_lag + y9_lag + y10_lag, data = data_t,
                                         mod_type = "triple", niter = niter, nburn = nburn, nthin = nthin,
                                        learn_a_xi = FALSE, learn_a_tau = FALSE,
                                        learn_c_xi = FALSE, learn_c_tau = FALSE)
                  
                  # Non-hierarchical triple gamma
                  res_TG <- shrinkTVP(y1 ~ y1_lag + y2_lag + y3_lag + y4_lag + y5_lag + y6_lag + y7_lag + y8_lag + y9_lag + y10_lag, data = data_t,
                                      mod_type = "triple", niter = niter, nburn = nburn, nthin = nthin,
                                      learn_kappa2_B = FALSE, learn_lambda2_B = FALSE,
                                      learn_a_xi = FALSE, learn_a_tau = FALSE,
                                      learn_c_xi = FALSE, learn_c_tau = FALSE)
                  
                  # Hierarchical horseshoe
                  res_H_HS <- shrinkTVP(y1 ~ y1_lag + y2_lag + y3_lag + y4_lag + y5_lag + y6_lag + y7_lag + y8_lag + y9_lag + y10_lag, data = data_t,
                                      mod_type = "triple", niter = niter, nburn = nburn, nthin = nthin,
                                      learn_a_xi = FALSE, learn_a_tau = FALSE,
                                      learn_c_xi = FALSE, learn_c_tau = FALSE,
                                      a_xi = 0.5, a_tau = 0.5, c_xi = 0.5, c_tau = 0.5)
                  
                  # Non-hierarchical horseshoe
                  res_HS <- shrinkTVP(y1 ~ y1_lag + y2_lag + y3_lag + y4_lag + y5_lag + y6_lag + y7_lag + y8_lag + y9_lag + y10_lag, data = data_t,
                                      mod_type = "triple", niter = niter, nburn = nburn, nthin = nthin,
                                      learn_kappa2_B = FALSE, learn_lambda2_B = FALSE,
                                      learn_a_xi = FALSE, learn_a_tau = FALSE,
                                      learn_c_xi = FALSE, learn_c_tau = FALSE,
                                      a_xi = 0.5, a_tau = 0.5, c_xi = 0.5, c_tau = 0.5)
                  
                  # Fully hierarchical double gamma
                  res_FH_DG <- shrinkTVP(y1 ~ y1_lag + y2_lag + y3_lag + y4_lag + y5_lag + y6_lag + y7_lag + y8_lag + y9_lag + y10_lag, data = data_t,
                                        niter = niter, nburn = nburn, nthin = nthin,
                                        hyperprior_param = list(nu_tau = 1, nu_xi = 1))
                  
                  # Hierarchical double gamma
                  res_H_DG <- shrinkTVP(y1 ~ y1_lag + y2_lag + y3_lag + y4_lag + y5_lag + y6_lag + y7_lag + y8_lag + y9_lag + y10_lag, data = data_t,
                                        niter = niter, nburn = nburn, nthin = nthin,
                                        learn_a_xi = FALSE, learn_a_tau = FALSE,
                                        a_xi = 0.1, a_tau = 0.1)
                  
                  # Non-hierarchical double gamma
                  res_DG <- shrinkTVP(y1 ~ y1_lag + y2_lag + y3_lag + y4_lag + y5_lag + y6_lag + y7_lag + y8_lag + y9_lag + y10_lag, data = data_t,
                                        niter = niter, nburn = nburn, nthin = nthin,
                                        learn_a_xi = FALSE, learn_a_tau = FALSE,
                                        a_xi = 0.1, a_tau = 0.1,
                                      learn_kappa2_B = FALSE, learn_lambda2_B = FALSE)
                  
                  # Hierarchical Lasso
                  res_H_LS <- shrinkTVP(y1 ~ y1_lag + y2_lag + y3_lag + y4_lag + y5_lag + y6_lag + y7_lag + y8_lag + y9_lag + y10_lag, data = data_t,
                                        niter = niter, nburn = nburn, nthin = nthin,
                                        learn_a_xi = FALSE, learn_a_tau = FALSE,
                                        a_xi = 1, a_tau = 1)
                  
                  # Non-hierarchical Lasso
                  res_LS <- shrinkTVP(y1 ~ y1_lag + y2_lag + y3_lag + y4_lag + y5_lag + y6_lag + y7_lag + y8_lag + y9_lag + y10_lag, data = data_t,
                                        niter = niter, nburn = nburn, nthin = nthin,
                                        learn_a_xi = FALSE, learn_a_tau = FALSE,
                                        a_xi = 1, a_tau = 1,
                                        learn_kappa2_B = FALSE, learn_lambda2_B = FALSE)
                  
                  # Ridge regression
                  res_FV <- shrinkTVP(y1 ~ y1_lag + y2_lag + y3_lag + y4_lag + y5_lag + y6_lag + y7_lag + y8_lag + y9_lag + y10_lag, data = data_t,
                                        mod_type = "ridge", niter = niter, nburn = nburn, nthin = nthin)
                  
                  lpds_res <- c(LPDS(res_FH_TG, data_test),
                                LPDS(res_H_TG, data_test),
                                LPDS(res_TG, data_test),
                                LPDS(res_H_HS, data_test),
                                LPDS(res_HS, data_test),
                                LPDS(res_FH_DG, data_test),
                                LPDS(res_H_DG, data_test),
                                LPDS(res_DG, data_test),
                                LPDS(res_H_LS, data_test),
                                LPDS(res_LS, data_test),
                                LPDS(res_FV, data_test))
                  
                  rm(list = ls()[!ls() %in% c("lpds_res", "us_data")])
                  
                  return(lpds_res)
                }
stopCluster(cl)
warnings()
cumu_lpds <- apply(lpds, 1, cumsum)
color <- c(rep("cyan3", 3),
           rep("firebrick3", 2),
           rep("forestgreen", 3),
           rep("yellow2",2),
           "black")
lty <- c(1:3, 1:2, 1:3, 1:2, 1)

# Plot results
par(mar=c(6,4,1,1))
colnames(cumu_lpds) <- c("Fully hierarchical NGG",
                         "Hierarchical NGG",
                         "NGG",
                         "Hierarchical Horseshoe",
                         "Horseshoe",
                         "Fully hierarchical NG",
                         "Hierarchical NG",
                         "NG",
                         "Hierarchical Lasso",
                         "Lasso",
                         "Ridge Regression")

matplot(cumu_lpds, type = "l", ylab = "Cumulative LPDS",
        xaxt = "n", xlab ="", col = color, lty = lty, lwd = 1.5)

# Extract labels from time series
labs = as.yearqtr(time(mym))[T0:Tmax + 1][c(FALSE, TRUE)]

# Create custom axis label
axis(1, at = (1:length(T0:Tmax))[c(FALSE,TRUE)], labels = FALSE)
text(x=(1:length(T0:Tmax))[c(FALSE,TRUE)],
     y=par()$usr[3]-0.05*(par()$usr[4]-par()$usr[3]),
     labels=labs, srt=45, adj=1, xpd=TRUE)
# Add legend
legend("topright", colnames(cumu_lpds), col= color,
       lty = lty, lwd = 1.5, bty = "n", cex = 0.8)


