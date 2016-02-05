## Set up
library(rstan)
library(arm)
library(loo)
library(abind)
library(assertthat)

rstan_options(auto_write = TRUE)

## needed on cluster
cores <- as.numeric(Sys.getenv("NSLOTS"))

if(is.na(cores))
    cores <- 2

options(mc.cores = cores)

set.seed(242345)

expose_stan_functions("logistic_impsampling_function.stan")

source("utils.R")

# Define the hyperparameters of the model
mu_a <- c(.5, -.2)
sigma_a <- c(.1, .1)
Sigma_a <- diag(sigma_a^2)
K <- length(mu_a)
b <- -.1

# Set up the experimental conditions
J <- 50
T <- 13
x <- (0:(T-1)) /T

# Simulate the data
a <- mvrnorm(J, mu_a, Sigma_a)
y <- array(NA, c(J,T))
n_trials <- 20
for (t in 1:T)
  y[,t] <- rbinom(J, n_trials, invlogit(a[,1] + a[,2]*x[t] + b*x[t]^2))

# Set up the conditions for the external data
J_prime <- 200
T_prime <- 13
x_prime <- (0:(T_prime-1)) / T_prime

# Simulate external data and compute averages
delta <- c(.1, .1)
a_prime <- mvrnorm(J_prime, mu_a + delta, Sigma_a)
y_prime <- array(NA, c(J_prime,T_prime))
n_trials_prime <- 20
for (t in 1:T_prime)
  y_prime[,t] <- rbinom(J_prime, n_trials_prime, invlogit(a_prime[,1] + a_prime[,2]*x_prime[t] + b*x_prime[t]^2))
y_prime_bar <- colMeans(y_prime)
true <- c(mu_a, b, sigma_a, delta)

# Set noninformative priors
K_phi <- 5
mu_phi_p <- rep(0,K_phi)
Sigma_phi_p <- diag(K_phi)
mu_delta_p <- rep(0,K)
Sigma_delta_p <- diag(K)

# Initialize the pseudo-prior for phi
mu_phi_g <- mu_phi_p
Sigma_phi_g <- Sigma_phi_p

parameters_ref <- c("mu_a", "b", "sigma_a", "delta")

# Fit the model to the complete data
fit_all <- stan("logistic_all.stan", data=nlist(T), iter=1000, chains=4)
print(fit_all, pars=parameters_ref)
sims_all <- extract(fit_all)
estimate_all <- colMeans(cbind(sims_all$mu_a, sims_all$b, sims_all$sigma_a, sims_all$delta))
sum_all <- summary(fit_all, pars=parameters_ref)

# Fit the model once just to local data
model_local <- stan("logistic_local.stan", data=nlist(T), chains=0)
fit_local <- stan(fit=model_local, data=nlist(T), warmup=2000, iter=2000+500, chains=4)

parameters_local <- exclude(parameters_ref, "delta")

sims <- extract(fit_local, c(parameters_local, "phi", "log_weight"))
print(fit_local, pars=c("mu_a", "b", "sigma_a", "lp__"))
##rstan::traceplot(fit_local, pars=c("mu_a", "b", "sigma_a"))
estimate_local <- c(colMeans(cbind(sims$mu_a, sims$b, sims$sigma_a)), NA, NA)

# Choose different starting points
n_rep <- 3
mu_delta_inits <- array(.5*rnorm(n_rep*K), c(n_rep,K))
Sigma_delta_g_init <- diag(K)

# set all EP algo parameters:
J_tilde <- 1000
n_phi_update <- 10
n_delta_update <- 10
n_save <- 100
chains <- 4
resample_steps <- 15
resample_lower <- 25
model_impsampling_rng <- logistic_impsampling_rng
n_warmup <- 2000
parameters <- rownames(sum_all$summary)

# run EP
source("extrap_ep_routine.R")

pdf_file <- paste0("extrap_logistic_", IS_weights,".pdf")

# Plot the parameter estimates and importance sampling efficiency as the algorithm runs
pdf(pdf_file, height=7, width=7)
par(oma=c(0,0,4,0), mfrow=c(3,3), mar=c(3,2,2,.5), mgp=c(1.5,.5,0), tck=-.02)
n_par <- dim(saved)[[3]]
for (k in 1:n_par){
  upper <- saved[,,k,1] + saved[,,k,2]
  lower <- saved[,,k,1] - saved[,,k,2]
  last_half <- (n_step/2 + 1):n_step
  plot(c(1, n_step), range(upper[,last_half], lower[,last_half], estimate_all[k], #estimate_correct[k],
                           estimate_local[k], na.rm=TRUE), type="n", ylab="", xlab=if (k>6) "Step of algorithm" else "", bty="l")
  for (i_rep in 1:n_rep)
    polygon(c(1:n_step, n_step:1), c(upper[i_rep,], rev(lower[i_rep,])), col="lightgray", border=NA)
  mtext(paste("Updating of", dimnames(saved)[[3]][k]), 3, 0, cex=.75)
  lines(c(1, n_step), rep(estimate_local[k], 2), col="red")
  lines(c(1, n_step), rep(estimate_all[k], 2), col="green")
  ## no correct estimate available
#  lines(c(1, n_step), rep(estimate_correct[k], 2), col="blue")
  for (i_rep in 1:n_rep)
    lines(1:n_step, saved[i_rep,,k,1])
}
plot(c(1, n_step), c(0,1.6), yaxs="i", type="n", ylab="", xlab="Step of algorithm", bty="l", main=expression(paste("Pareto shape ", hat(k), " of importance weights")), cex.main=1)
for (i_rep in 1:n_rep) lines(1:n_step, pareto[i_rep,], lwd=.5)
abline(1,0,lty=3)
abline(0.5,0,lty=3)
plot(c(1, n_step), range(0, imp_efficiency[i_rep,])*1.05, yaxs="i", type="n", ylab="", xlab="Step of algorithm", bty="l", main="Efficiency of importance samples", cex.main=1, font.main=1)
for (i_rep in 1:n_rep) lines(1:n_step, imp_efficiency[i_rep,], lwd=.5)
mtext(paste("Hierarchical logistic example:  Posterior mean +/- sd from EP algorithm from", n_rep, "starting points\n(Red lines show estimate from local data, green uses complete data)"), line=1, side=3, outer=TRUE, cex=.8)
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


save(list=ls(), file=paste0("extrap_logistic_", IS_weights, ".rda"))

sessionInfo()

