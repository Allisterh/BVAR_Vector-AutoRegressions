library(coda)
install.packages("rstan")
library(rstan)
options(mc.cores = parallel::detectCores())

chains=8
it_burnin = 10000
it_sample = 30000

sm = stan_model('Desktop/Project/MS_VAR_wishart.stan',auto_write=TRUE)
#sm = stan_model('Stan/MS_VAR_lkj.stan',auto_write=TRUE)

######################################################################################

n=nrow(dsg)

######################################################################################

n_dim=ncol(dsg)
n_AR=1
n_reg=2
reg_m=TRUE
reg_s=TRUE

stanfit = sampling(sm,data=list(T=n,dim=n_dim,y=dsg,ARdim=n_AR,nreg=n_reg,mean_reg=reg_m,sigma_reg=reg_s,Q_alpha=as.array(rep(1,n_reg)),
                                mu_mean=0,mu_sd=0.2,phi_mean=0,phi_sd=1,eye=diag(n_dim),L=chol(diag(n_dim)),nu=n_dim+1),
                   iter=(it_burnin+it_sample), warmup = it_burnin, chains=chains,
                   pars=c("Q","mu","phi","sigma","S"))

# stanfit = sampling(sm,data=list(T=n,dim=n_dim,y=Y,ARdim=n_AR,nreg=n_reg,mean_reg=reg_m,sigma_reg=reg_s,Q_alpha=as.array(rep(1,n_reg)),
#                    mu_mean=0,mu_sd=0.2,phi_mean=0,phi_sd=1,eta=1,gamma_alpha=2,gamma_beta=0.1),
#                    iter=(it_burnin+it_sample), warmup = it_burnin, chains=chains,
#                    pars=c("Q","mu","phi","sigma","S"))

######################################################################################
# Label-switching

pars =mcmc.list(lapply(1:ncol(stanfit), function(x) mcmc(as.array(stanfit)[,x,])))

source("Desktop/Project/labsw_2reg.R")
pars_ls = labsw_2reg(pars,chains,n_dim,n_AR,n_reg,reg_m,reg_s)

######################################################################################
# Results

options(scipen=999)
summary(pars_ls$pars)$statistics[1:pars_ls$index[1],]
plot(pars_ls$pars[,c(1,2,pars_ls$index[2],pars_ls$index[2]+1)])

states=as.matrix(pars_ls$pars)[,(pars_ls$index[1]+1):pars_ls$index[3]]
plot(colMeans(states),ylab="Mean state")
plot(round(colMeans(states),0),ylab="Most likely state")
ggplot(data, aes(x,-y-axis) ) +geom_point() + geom_smooth(method = "lm")
