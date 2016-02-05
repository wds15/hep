
## numerically stable summation of logs on the natural scale
log_sum_exp <- function(x){ 
   xmax <- which.max(x) 
   log1p(sum(exp(x[-xmax]-x[xmax])))+x[xmax] 
} 

## log truncation of importance weights
log_trim <- function(lw) {
    S <- length(lw)
    m <- log_sum_exp(lw) - log(S)
    pmin(lw, 0.5 * log(S) + m)
}

## scale weights on average to 1
scale_log_weights <- function(lw) {
    wmean <- log_sum_exp(lw)
    lw - wmean
}

## exclude elements from a list
exclude <- function(set, elems, keys=set) {
    set[-match(elems, keys)]
}

## numerically stable antilog of weights
exp_norm <- function(log_w) {
  return(exp(log_w - max(log_w)))
}

## column wise variance
colVars <- function(a) {
    n <- dim(a)[[1]]
    c <- dim(a)[[2]]
    return(.colMeans(((a - matrix(.colMeans(a, n, c), nrow = n, ncol = c, byrow = TRUE)) ^ 2), n, c) * n / (n - 1))
}


## take a posterior and resample by taking a fraction of Seff samples
## without/with replacement; if any if the log_weights underflows, the
## resampling is done one by one
resample_post <- function(post, log_weights, frac=0.5, lower=min(25, length(log_weights)), replace=FALSE) {
    weights <- exp_norm(log_weights)
    assert_that(NROW(post[[1]]) == length(weights))
    S <- length(weights)
    Seff <- efficiency_weights(log_weights, TRUE) * S
    Nresamp <- max(lower, round(frac*Seff))
    cat("Sampling", Nresamp, "draws ( Seff =", round(Seff), ")\n")
    ## sample indices without replacement
    if(!replace & any(weights==0)) {
        ## do it one-by-one for numerical stability reasons
        cand_set <- 1:S
        ind <- c()
        while(length(ind) != Nresamp) {
            ind_draw <- sample.int(length(cand_set), 1, prob=exp_norm(log_weights))
            ind <- c(ind, cand_set[ind_draw])
            log_weights <- log_weights[-ind_draw]
            cand_set <- cand_set[-ind_draw]
        }
    } else {
        ind <- sample.int(S, Nresamp, prob=weights, replace=replace)
    }
    lapply(post, asub, idx=list(ind), dim=1)
}

## calculate approximatley the fraction of independent samples the
## weights represent
efficiency_weights <- function(w, log=FALSE) {
    S <- length(w)
    if(log) {
        w <- w - log_sum_exp(w) + log(S)
        efficiency <- exp(log(S) - log_sum_exp(2*w))
    } else {
        w_norm <- w/mean(w)
        efficiency <- length(w)/sum(w_norm^2)
    }
    return(efficiency)
}

## EP utility function
ep_diff <- function(Sigma_inv_1, Sigma_inv_2, Sigma_inv_mu_1, Sigma_inv_mu_2) {
  pd_check <- function(a) return(all(eigen(a)$values>0))
  pd <- FALSE
  n <- -1
  while(!pd){
    n <- n + 1
    Sigma_inv <- Sigma_inv_2 - Sigma_inv_1/2^n
    pd <- pd_check(Sigma_inv)
  }
  Sigma_inv_mu <- Sigma_inv_mu_2 - Sigma_inv_mu_1/2^n
  return(nlist(Sigma_inv, Sigma_inv_mu, n))
}

## extract from a posterior given as as list a specifc draw. Assumes
## that the first dimension of each list entry is the iteration.
extract_draw <- function(sims, draw) lapply(sims, asub, idx=draw, dim=1)

## applies a function over each entry of the posterior if
## vectorized=FALSE; for vectorized=TRUE the function is assumed to
## perform the simulation in a single sweep. Note that all arguments
## to the function are automatically deduced from it's formals and
## that all arguments which are not in the sims list are searched in
## the global environment.
sim_posterior <- function(sims, fun, vectorized=FALSE) {
    args <- setdiff(names(formals(fun)), "seed")

    from_draw <- intersect(args, names(sims))
    from_env  <- setdiff(args, names(sims))

    sims <- sims[from_draw]
    aux <- mget(from_env, envir=parent.frame())

    if(!vectorized) {
        S <- NROW(sims[[1]])
        calc_draw <- function(i) do.call(fun, c(aux, extract_draw(sims, i)))
        res_type <- calc_draw(1)
        return(t(vapply(1:S, calc_draw, res_type)))
    } else {
        return(do.call(fun, c(sims, aux)))
    }
}

## calculates reweighted mean and variance given log-weights
weighted_mean_and_var <- function (a, lw) {
    N <- length(lw)
    ess <- N*efficiency_weights(lw, log=TRUE)
    if(ess <= 1) {
        ## for values <= 1 the ess correction of the variance below
        ## collapses, hence correct for this
        warning(paste("ESS =", ess, "artifically increased to", 1.1))
        ess <- 1.1
    }
    w <- exp_norm(lw)
    if (is.matrix(a)){
        K <- ncol(a)
        a_mean <- rep(NA, K)
        a_var  <- array(NA, c(K,K))
        for (k in 1:K)
            a_mean[k] <- weighted.mean(a[,k], w)
        for (k1 in 1:K)
            for (k2 in 1:K)
                a_var[k1,k2] <- weighted.mean((a[,k1] - a_mean[k1])*(a[,k2] - a_mean[k2]), w)
        a_var <- a_var*ess/(ess - 1)
        a_sd  <- sqrt(diag(a_var))
    } else {
        a_mean <- weighted.mean(a, w)
        a_var  <- weighted.mean((a - a_mean)^2, w)
        a_var  <- a_var*ess/(ess - 1)
        a_sd   <- sqrt(a_var)
    }
    return(list(mean=a_mean, var=a_var, sd=a_sd))
}

