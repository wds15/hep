
# Set up the EP procedure
iter <- round(100*2^(((1:n_phi_update)-1)/2))
n_step <- n_phi_update*n_delta_update
saved <- array(NA, c(n_rep, n_step, length(parameters), 2))
dimnames(saved) <- list(1:n_rep, 1:n_step, parameters, c("mean","sd"))
imp_efficiency <- array(NA, c(n_rep, n_step))
pareto <- array(NA, c(n_rep, n_step))

# Saving things to do R-hat
resampled_sims <- array(NA, c(n_step*n_save, n_rep, length(parameters)))

iterMult <- ceiling(chains / cores)

verbose <- TRUE

n_max <- round(2*2^(((1:n_phi_update)-1)/2))

## importance reweighting method can be
## IS with vanllia importance weights
## TIS with truncated importance weights
## PSIS with Pareto smoothed importance weights
##IS_weights <- "TIS"
IS_weights <- "PSIS"
##IS_weights <- "IS"

# Our algorithm!
for (i_rep in 1:n_rep){
    cat("Replication", i_rep, "of", n_rep, "\n")
    mu_phi_g <- mu_phi_p
    Sigma_phi_g <- Sigma_phi_p
    mu_delta_g <- mu_delta_inits[i_rep,]
    Sigma_delta_g <- Sigma_delta_g_init
    step <- 0
    for (loop in 1:n_phi_update){
        cat("  Loop", loop, "of", n_phi_update, "\n")
        # Fit data y
        fit_local <- stan(fit=model_local, data=nlist(T), cores=cores, warmup=n_warmup, iter=n_warmup + iterMult*iter[loop], chains=chains)
        n_sims <- chains * iterMult*iter[loop]
        local_sims <- extract(fit_local, c(parameters_local, "phi", "log_weight"))
        local_sims_parameters <- extract(fit_local, parameters_local, inc_warmup=FALSE, permuted=FALSE)
        local_sims_parameters <- matrix(local_sims_parameters, nrow=n_sims, ncol=dim(local_sims_parameters)[3])
        
        # Throw in external data y_prime_bar
        for (i in 1:n_delta_update) {
            step <- step + 1
            cat("Replication", i_rep, "loop", loop, "delta update", i, "of", n_delta_update, "\n")
            gen <- sim_posterior(local_sims, model_impsampling_rng)
            generated <- list(delta=gen[,1:K,drop=FALSE], log_weight=gen[,K+1])
            log_weight <- local_sims$log_weight + generated$log_weight
            weight <- exp_norm(log_weight)
            psis <- psislw(log_weight)
            pareto_k <- psis$pareto_k
            if(IS_weights == "IS") {
                ## vanilla importance sampling weights
                log_tweight <- log_weight
            } else if(IS_weights == "TIS") {
                ## truncated importance sampling weights
                log_tweight <- log_trim(log_weight)
            } else {
                ## by default, use Pareto smoothed weights
                log_tweight <- psis$lw
                if(any(!is.finite(log_tweight))) {
                    warning("Replication ", i_rep, " loop", loop, " delta update ", i, " of ", n_delta_update, " switch to log_trim.")
                    cat("Replication", i_rep, "loop", loop, "delta update", i, "of", n_delta_update, "switch to log_trim.\n")
                    log_tweight <- log_trim(log_weight)
                }
            }
            tweight <- exp_norm(log_tweight)
            eff  <- efficiency_weights(log_weight , TRUE)
            imp_efficiency[i_rep,step] <- eff
            pareto[i_rep,step] <- pareto_k
            
            cat("Seff     =", round(eff *n_sims, 2), "/", round(100*eff, 2), "%\npareto_k =", pareto_k, "\n")
      
            draws <- sample.int(n_sims, n_save, replace=TRUE, prob=tweight)
            index <- (step-1)*n_save + (1:n_save)
            resampled_sims[index,i_rep,] <- cbind(local_sims_parameters[draws,], generated$delta[draws,])

            ## Do resampling without replacement during the first
            ## steps (defined by resample_steps) as in the beginning
            ## weights can be unstable as theses collapse on very few
            ## elements and hence leading to an underestimated
            ## variance.  Resampling tends to over-estimate variance,
            ## but this is better than crashing.
            local_sims_cand <- c(local_sims, list(delta=generated$delta))
            if(step <= resample_steps) {
                ## resampling without replacement is done with the
                ## "vanilla" weights
                local_sims_cand <- resample_post(local_sims_cand, log_weight, lower=resample_lower)
                ## now things are equally weighted
                log_tweight <- rep(0, NROW(local_sims_cand[[1]]))
            }
            
            post_est <- lapply(local_sims_cand[parameters_ref], weighted_mean_and_var, lw=log_tweight)

            Sigma_delta_g_cand <- matrix(post_est$delta$var, K, K)

            if(det(Sigma_delta_g_cand) <= 0) {
                n_delta_g <- diag(Sigma_delta_p) / diag(Sigma_delta_g_cand)
                scale_delta_g <- pmax(n_delta_g / n_max[loop], 1)
                Sigma_delta_g_cand <- Sigma_delta_g_cand + diag(diag(Sigma_delta_g_cand) * (scale_delta_g-1), K, K)
            }
      
            post_mean <- unlist(lapply(post_est, "[[", "mean"))
            post_sd   <- unlist(lapply(post_est, "[[", "sd"))

            saved[i_rep,step,,1] <- post_mean
            saved[i_rep,step,,2] <- post_sd
            mu_delta_g <- post_est$delta$mean
            Sigma_delta_g <- Sigma_delta_g_cand
        }

        if(verbose) {
            cat("Current estimates:\n")
            print(saved[i_rep,step,,1])
            cat("Current sds:\n")
            print(saved[i_rep,step,,2])
            cat("BEFORE EP update:\n")
            cat("mu_phi_g =", mu_phi_g,"\n")
            cat("sqrt(diag(Sigma_phi_g)) =", sqrt(diag(Sigma_phi_g)),"\n")
        }
        
        # EP calculations
        zeros <- rep(0, n_sims)
        temp1 <- weighted_mean_and_var(local_sims$phi, zeros)
        Sigma_inv_1<- solve(temp1$var)
        Sigma_inv_mu_1<- solve(temp1$var) %*% temp1$mean
        
        temp2 <- weighted_mean_and_var(local_sims_cand$phi, log_tweight)
        Sigma_inv_2 <- solve(temp2$var) + solve(Sigma_phi_g)
        Sigma_inv_mu_2 <- solve(temp2$var) %*% temp2$mean + solve(Sigma_phi_g) %*% mu_phi_g
            
        diff <- ep_diff(Sigma_inv_1, Sigma_inv_2, Sigma_inv_mu_1, Sigma_inv_mu_2)
        mu_phi_g <- as.vector(solve(diff$Sigma_inv, diff$Sigma_inv_mu))
        Sigma_phi_g <- solve(diff$Sigma_inv)
        
        ## calculate roughly how many observations the proposal is
        ## worth in units of the prior reference scale
        n_phi_g <- diag(Sigma_phi_p) / diag(Sigma_phi_g)
        scale_var_g <- pmax(n_phi_g / n_max[loop], 1)
        Sigma_phi_g <- Sigma_phi_g + diag(diag(Sigma_phi_g) * (scale_var_g-1))
        
        if(verbose) {
            cat("AFTER EP update:\n")
            cat("mu_phi_g =", mu_phi_g,"\n")
            cat("sqrt(diag(Sigma_phi_g)) =", sqrt(diag(Sigma_phi_g)),"\n")
        }
    }
}

