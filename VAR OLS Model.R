# Estimating a VAR model in R
library(vars)
library(tseries)
library(tidyverse)
library(stargazer)
library(urca)
library(sandwich)
library(readr)
library(xts)
library(zoo)
library(forecast)
library(ggplot2)
library(fpp2)
library(quantmod)

# Need to check for integration properties of my series, check I(0) meaning they are stationary
# Need to choose optimal lag length, too many lags might lead to dimensionality problems with too many parameters to estimate. Too little, might lead to serial/auto correlation.
# We use Information Crtieria to tell us optimal lag length
# Package Vars esimates parameters equation by equation using OLS
# Need to make sure model is stable by checking eigenvalues being less than 1, if it is model is stationary
# Granger causality test
# Generate IRFs and Variance decomposition

head(Complete_Quarterly_Dataset)
tail(Complete_Quarterly_Dataset)
str(Complete_Quarterly_Dataset)
#declaring the data is time series
Date = Complete_Quarterly_Dataset$Date
y1<- ts(Complete_Quarterly_Dataset$y1, start=c(1960,2), frequency = 4)
y2<- ts(Complete_Quarterly_Dataset$y2, start=c(1960,2), frequency = 4)
y3<- ts(Complete_Quarterly_Dataset$y3, start=c(1960,2), frequency = 4)
y4<- ts(Complete_Quarterly_Dataset$y4, start=c(1960,2), frequency = 4)
y5<- ts(Complete_Quarterly_Dataset$y5, start=c(1960,2), frequency = 4)
y6<- ts(Complete_Quarterly_Dataset$y6, start=c(1960,2), frequency = 4)
y7<- ts(Complete_Quarterly_Dataset$y7, start=c(1960,2), frequency = 4)
y8<- ts(Complete_Quarterly_Dataset$y8, start=c(1960,2), frequency = 4)
y9<- ts(Complete_Quarterly_Dataset$y9, start=c(1960,2), frequency = 4)
y10<- ts(Complete_Quarterly_Dataset$y10, start=c(1960,2), frequency = 4)
plot.ts(y1)

# Checking I(0) using Augmented Dickey-Fuller Test (ADF)
adfy1 <- adf.test(y1)
print(adfy1)
# Data is already stationary as I applied 400*Ln(y_t/y_t-1) in Excel.
# I will now go variable by variable checking and correcting non-stationarity.
adfy2 <- adf.test(Complete_Quarterly_Dataset$y2)
print(adfy2)
plot.ts(y2)
adfy3 <- adf.test(Complete_Quarterly_Dataset$y3)
print(adfy3)
adfy4 <- adf.test(Complete_Quarterly_Dataset$y4)
print(adfy4)
adfy5 <- adf.test(Complete_Quarterly_Dataset$y5)
print(adfy5)
plot.ts(y5)
adfy6 <- adf.test(Complete_Quarterly_Dataset$y6)
print(adfy6)
adfy7 <- adf.test(Complete_Quarterly_Dataset$y7)
print(adfy7)
adfy8 <- adf.test(Complete_Quarterly_Dataset$y8)
print(adfy8)
adfy9 <- adf.test(Final_Quarterly_Dataset$y9)
print(adfy9)
adfy10 <- adf.test(Final_Quarterly_Dataset$y10)
print(adfy10)
# Choosing optimal lag length based on information criteria. Too few lags you will have correlation which means information is still included in the error term that hasn't been extracted yet. 
mydata = cbind(y1, y2, y3, y4, y5, y6, y7, y8, y9, y10)
colnames(mydata) <- c("GDP", "PFCE", "IPI", "UR", "M3", "FFER", "CPI", "CI", "LTIR","CR")
lag <- VARselect(mydata)
lag
# Monte Carlo evidence suggests that the BIC or SC penalty function works well in designs calibrated to macroeconomic data (Stock and Watson, 2011) 
var1 <- VAR(mydata, p=1, type = "none")
summary(var1)
pairs(Complete_Quarterly_Dataset[,2:11])
install.packages("stargazer")
# Stargazer is better than the summary function as you can read equations vertically.
stargazer(var1[["varresult"]], type='text')

# Checking if the model is stable by checking eigenvalues are less than 1.
roots(var1, modulus = TRUE)
# Eigenvalues are less than 1 so the model is stable, this means that forecasts will converge towards the unconditional mean of the variables in the model, and the MSE matrix of the forecast errors will converge towards the unconditional variance of the variables in the model.

# Granger causality test
grangertest(y4, y1, order = 1)
grangery10 <- causality(var1, cause = "y10")
grangery10
# We have a large p-value, more than 5%, so we can't reject the null hypothesis.
grangery8 <- causality(var1,cause = "y8")
grangery8
# Small p-value so we can't reject the null hypothesis

#IRFs 
irf1 <- irf(var1, impulse = "y8", response = "y7", n.ahead = 20, boot = TRUE, run=200, ci=0.95)
plot(irf1, ylab="Inflation", main = "Inflation response to Confidence Indicator shock")
summary(var1)$varresult$y8$sigma
irf2 <- irf(var1, impulse = "y6", response = "y4", n.ahead = 20, boot = TRUE, run=200, ci=0.95)
plot(irf2, ylab="Unemployment", main = "Unemployment response to Fed Funds Effective Rate shock")
summary(var1)$varresult$y6$sigma
irf3 <- irf(var1, impulse = "y9", response = "y4", n.ahead = 20, boot = TRUE, run=200, ci=0.95)
plot(irf3, ylab="Unemployment", main = "Unemployment response to 10-Year Bond Yield shock")
summary(var1)$varresult$y9$sigma
# Variance decomposition, changing the structure a bit for visibility
quartz(width=12,height=13)
plot(fevd(var1, n.ahead=10))
fevd1 <- fevd(var1, n.ahead=10)
plot.varfevd  <-function (x, plot.type = c("multiple", "single"), names = NULL,
                          main = NULL, col = NULL, ylim = NULL, ylab = NULL, xlab = NULL,
                          legend = NULL, names.arg = NULL, nc, mar = par("mar"), oma = par("oma"),
                          addbars = 1, ...)
{
  K <- length(x)
  ynames <- names(x)
  plot.type <- match.arg(plot.type)
  if (is.null(names)) {
    names <- ynames
  }
  else {
    names <- as.character(names)
    if (!(all(names %in% ynames))) {
      warning("\nInvalid variable name(s) supplied, using first variable.\n")
      names <- ynames[1]
    }
  }
  nv <- length(names)
  #    op <- par(no.readonly = TRUE)
  ifelse(is.null(main), main <- paste("FEVD for", names), main <- rep(main,
                                                                      nv)[1:nv])
  ifelse(is.null(col), col <- gray.colors(K), col <- rep(col,
                                                         K)[1:K])
  ifelse(is.null(ylab), ylab <- rep("Percentage", nv), ylab <- rep(ylab,
                                                                   nv)[1:nv])
  ifelse(is.null(xlab), xlab <- rep("Horizon", nv), xlab <- rep(xlab,
                                                                nv)[1:nv])
  ifelse(is.null(ylim), ylim <- c(0, 1), ylim <- ylim)
  ifelse(is.null(legend), legend <- ynames, legend <- legend)
  if (is.null(names.arg))
    names.arg <- c(paste(1:nrow(x[[1]])), rep(NA, addbars))
  plotfevd <- function(x, main, col, ylab, xlab, names.arg,
                       ylim, ...) {
    addbars <- as.integer(addbars)
    if (addbars > 0) {
      hmat <- matrix(0, nrow = K, ncol = addbars)
      xvalue <- cbind(t(x), hmat)
      barplot(xvalue, main = main, col = col, ylab = ylab,
              xlab = xlab, names.arg = names.arg, ylim = ylim,
              legend.text = legend, ...)
      abline(h = 0)
    }
    else {
      xvalue <- t(x)
      barplot(xvalue, main = main, col = col, ylab = ylab,
              xlab = xlab, names.arg = names.arg, ylim = ylim,
              ...)
      abline(h = 0)
    }
  }
  if (plot.type == "single") {
    #        par(mar = mar, oma = oma)
    #        if (nv > 1)
    #            par(ask = TRUE)
    for (i in 1:nv) {
      plotfevd(x = x[[names[i]]], main = main[i], col = col,
               ylab = ylab[i], xlab = xlab[i], names.arg = names.arg,
               ylim = ylim, ...)
    }
  }
  else if (plot.type == "multiple") {
    if (missing(nc)) {
      nc <- ifelse(nv > 4, 2, 1)
    }
    nr <- ceiling(nv/nc)
    par(mfcol = c(nr, nc), mar = mar, oma = oma)
    for (i in 1:nv) {
      plotfevd(x = x[[names[i]]], main = main[i], col = col,
               ylab = ylab[i], xlab = xlab[i], names.arg = names.arg,
               ylim = ylim, ...)
    }
  }
  #    on.exit(par(op))
}
quartz(width=15,height=14)
layout(matrix(1:10, ncol=2))
plot.varfevd(fevd(var1, n.ahead = 10), plot.type = "single", col = 1:10)
# Breakdown of where the shock comes from, use legend on the side to interpret the graphs.
# Making predictions below, first column is the forecast
betapredic=predict(var1, n.ahead=4, ci=0.95, dumvar=NULL)
betapredic
fanchart(betapredic, names = "y1")
fanchart(betapredic, names = "y3")
fanchart(betapredic, names = "y4")
