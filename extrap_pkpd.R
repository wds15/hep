library(assertthat)
library(rstan)
library(reshape2)
library(ggplot2)
library(plyr)
library(parallel)
library(functional)
library(StanHeaders)
library(arm)
library(loo)
library(Rcpp)
library(abind)
theme_set(theme_bw())
source("utils.R")

rstan_options(auto_write = TRUE)

set.seed(24234)

## choose whatever is appropiate for your computing environment
## useful on laptop
##cores <- detectCores()
## needed on cluster
cores <- as.numeric(Sys.getenv("NSLOTS"))

options(mc.cores = cores)

expose_stan_functions("pkpd.stan")

# define true parameters
lalpha_0 <- log(50.)
  
lalpha_s <- log(42.)
lkappa   <- log(42.) - 2*log(10)
lEmax_s  <- log(0.4)
    
sigma_lalpha_0 <- 0.10
sigma_lalpha_s <- 0.15

sigma_y <- 0.2

J <- 100 ## from the first group the first half is placebo, the rest is on treatment 1
J_prime <- 50

## let's say we cover a year and measure monthly
x <- seq(0, 52, length=13)
T <- 13

## for simplicity assume the prime data set is at the same time-points
x_prime <- seq(0, 52, length=13)
T_prime <- 13

## number of simulations per draw of the posterior to get
## approximated log-lik weight
J_tilde <- 1e3

## number of simulations per draw of the posterior
J_post <- 1e3

## define weakly-informative prior
# phi is ordered as:
##  phi[1] <- lalpha_0;
##  phi[2] <- lalpha_s;
##  phi[3] <- lkappa;
##  phi[4] <- lEmax_s;
##  phi[5] <- log(sigma_lbva);
##  phi[6] <- log(sigma_lalpha_s);
##  phi[7] <- log(sigma_y);
mu_phi_p <- c(log(50), log(50), log(50) - 2, log(0.1), log(0.1), log(0.1), log(1))
K_phi <- length(mu_phi_p)
Sigma_phi_p <- diag(rep(5^2, K_phi))

## initialize proposal density with prior
mu_phi_g <- mu_phi_p
Sigma_phi_g <- Sigma_phi_p

K <- 1
delta <- array(0.2, dim=1)

mu_delta_p <- array(0, dim=1)
Sigma_delta_p <- matrix(5^2)

## simulated asymptotic values for each group
exp(lalpha_s) ## placebo t=inf
exp(lalpha_s + log(1 + exp(lEmax_s))) ## drug=1 t=inf
exp(lalpha_s + log(1 + exp(lEmax_s + delta))) ## drug=2 t=inf


## simulate patient specific parameters
J_tot <- J + J_prime
subj_tot <- data.frame(id=1:J_tot
                      ,prime=rep(c(0,1), times=c(J, J_prime))
                      ,DRUG=0
                      ,eta_lalpha_0=rnorm(J_tot, 0, 1)
                      ,eta_lalpha_s =rnorm(J_tot, 0, 1)
                   )

## to reduce simulation noise influence, we center the random effects
## exactly to 0
subj_tot <-  transform(subj_tot
                      ,eta_lalpha_0=eta_lalpha_0-mean(eta_lalpha_0)
                      ,eta_lalpha_s=eta_lalpha_s-mean(eta_lalpha_s))

subj_tot$DRUG[(J/2):J] <- 1
subj_tot$DRUG[(J+1):J_tot] <- 2

subj       <- subset(subj_tot, prime==0)
subj_prime <- subset(subj_tot, prime==1)

lysim <- sim_posterior(subj, evaluate_model, TRUE)
lysim_prime <- sim_posterior(subj_prime, evaluate_model, TRUE)

dimnames(lysim) <- list(id=1:J, T=1:T)
dimnames(lysim_prime) <- list(id=(1+J):J_tot, T=1:T_prime)

M <- arrange(melt(lysim), id, T)
Mp <- arrange(melt(lysim_prime), id, T)
M$x <- x[M$T]
Mp$x<- x_prime[Mp$T]

M <- merge(rbind(M,Mp), subj_tot, by="id")

## look at the simulated means
plm <- ggplot(M, aes(x, exp(value), id=id, colour=factor(DRUG))) +
    geom_line(alpha=0.3) + 
        stat_summary(aes(id=DRUG), fun.data = "mean_cl_boot", position=position_dodge(width=0.5))


## create noise which we recycle for all of the different delta
## simulations
ly_err       <- matrix(rnorm(J*T,             0, sigma_y), J,       T)
ly_prime_err <- matrix(rnorm(J_prime*T_prime, 0, sigma_y), J_prime, T_prime)
ly       <- lysim       + ly_err

ly_prime <- lysim_prime + ly_prime_err

ly_prime_bar <- colMeans(ly_prime)

y <- exp(ly)
y_prime <- exp(ly_prime)
y_prime_bar <- exp(ly_prime_bar)

## perform first fit & calculate importance weights along

model_local <- stan("pkpd.stan", data=nlist(T), chains=0)

DRUG <- subj_tot$DRUG
DRUG_prime <- 2

stan_fit_all <- stan(fit=model_local, data=list(fit_all=1, T=T),
                     chains=4, iter=500+500, warmup=500,
                     seed = sample.int(.Machine$integer.max, 1) )


parameters_ref <- c("lalpha_0", "lalpha_s", "lkappa", "lEmax_s", "sigma_lalpha_0", "sigma_lalpha_s", "sigma_y", "delta")
parameters_local <- exclude(parameters_ref, "delta")

print(stan_fit_all, pars=parameters_ref )
sims_all <- extract(stan_fit_all)
sum_all <- summary(stan_fit_all, pars=parameters_ref)
estimate_all <- sum_all$summary[,"mean"]


# Fit the model once just to local data

## make sure the following fits do not fit all data
fit_all <- 0
fit_local <- stan(fit=model_local, data=nlist(T),
                chains=4, iter=500+500, warmup=500,
                seed = sample.int(.Machine$integer.max, 1) )

sims <- extract(fit_local, c(parameters_local, "phi", "log_weight"))
print(fit_local, pars=parameters_local)
sum_local <- summary(fit_local, pars=parameters_local)
estimate_local <- c(sum_local$summary[,"mean"], NA)


# set all EP algo parameters:
n_rep <- 3
J_tilde <- 1000
n_phi_update <- 10
n_delta_update <- 10
n_save <- 100
chains <- 4
resample_steps <- 15
resample_lower <- 25
model_impsampling_rng <- pkpd_impsampling_rng
n_warmup <- 500
parameters <- rownames(sum_all$summary)

# Choose different starting points
mu_delta_inits <- array(.5*rnorm(n_rep*K), c(n_rep,K))
Sigma_delta_g_init <- diag(K)

# run EP
source("extrap_ep_routine.R")

pdf_file <- paste0("extrap_pkpd_", IS_weights,".pdf")

# Plot the parameter estimates and importance sampling efficiency as the algorithm runs
pdf(pdf_file, height=7, width=7)
par(oma=c(0,0,4,0), mfrow=c(3,3), mar=c(3,2,2,.5), mgp=c(1.5,.5,0), tck=-.02)
n_par <- dim(saved)[[3]]
dnames <- dimnames(saved)
for (k in 1:n_par){
    ## skip over sigma_y
    if(dnames[[3]][k] == "sigma_y")
        next;
    upper <- saved[,,k,1] + saved[,,k,2]
    lower <- saved[,,k,1] - saved[,,k,2]
    last_half <- (n_step/2 + 1):n_step
    plot(c(1, n_step), range(upper[,last_half], lower[,last_half], estimate_all[k],
                             estimate_local[k], na.rm=TRUE), type="n", ylab="", xlab=if (k>6) "Step of algorithm" else "", bty="l")
    for (i_rep in 1:n_rep)
        polygon(c(1:n_step, n_step:1), c(upper[i_rep,], rev(lower[i_rep,])), col="lightgray", border=NA)
    mtext(paste("Updating of", dimnames(saved)[[3]][k]), 3, 0, cex=.75)
    lines(c(1, n_step), rep(estimate_local[k], 2), col="red")
    lines(c(1, n_step), rep(estimate_all[k], 2), col="green")
    for (i_rep in 1:n_rep)
        lines(1:n_step, saved[i_rep,,k,1])
}
plot(c(1, n_step), c(0,1.6), yaxs="i", type="n", ylab="", xlab="Step of algorithm", bty="l", main=expression(paste("Pareto shape ", hat(k), " of importance weights")), cex.main=1)
for (i_rep in 1:n_rep) lines(1:n_step, pareto[i_rep,], lwd=.5)
abline(1,0,lty=3)
abline(0.5,0,lty=3)
plot(c(1, n_step), range(0, imp_efficiency[i_rep,])*1.05, yaxs="i", type="n", ylab="", xlab="Step of algorithm", bty="l", main="Efficiency of importance samples", cex.main=1, font.main=1)
for (i_rep in 1:n_rep) lines(1:n_step, imp_efficiency[i_rep,], lwd=.5)
mtext(paste("Hierarchical PK/PD example:  Posterior mean +/- sd from EP algorithm from", n_rep, "starting points\n(Red lines show estimate from local data, green uses complete data)"), line=1, side=3, outer=TRUE, cex=.8)
dev.off()

## convert PDF to PS
system(paste("pdftops", pdf_file))

# R-hat using the simulations of theta

monitor(resampled_sims, digits=2)

# R-hat using the analytic formula

saved_perm <- array(NA, dim(saved)[c(2,1,3)])
for (k in 1:n_par){
  saved_perm[,,k] <- t(saved[,,k,1])
}
monitor(saved_perm, digits=2)

save(list=ls(), file=paste0("extrap_pkpd_", IS_weights, ".rda"))

sessionInfo()

