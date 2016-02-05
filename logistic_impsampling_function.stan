functions {
  vector logistic_impsampling_rng(int J_prime, int T_prime, int K, vector y_prime_bar, vector x_prime, vector mu_delta_p, matrix Sigma_delta_p, vector mu_delta_g, matrix Sigma_delta_g, int J_tilde, vector mu_a, real b, vector sigma_a, int n_trials_prime) {
    vector[K+1] output;
    vector[K] delta;
    vector[T_prime] M_tilde;
    matrix[T_prime,T_prime] Sigma_tilde;
    real log_weight;
    delta <- multi_normal_rng(mu_delta_g, Sigma_delta_g);
    {
      vector[K] mu_a_tilde;
      matrix[J_tilde,K] a_tilde;
      int y_tilde[J_tilde,T_prime];
      matrix[J_tilde,T_prime] y_tilde_real;
      matrix[J_tilde,T_prime] y_tilde_centered;
      mu_a_tilde <- mu_a + delta;
      for (j in 1:J_tilde){
        for (k in 1:K)
          a_tilde[j,k] <- normal_rng(mu_a_tilde[k], sigma_a[k]);
        for (t in 1:T_prime){
          y_tilde[j,t] <- binomial_rng(n_trials_prime, inv_logit(a_tilde[j,1] + a_tilde[j,2]*x_prime[t] + b*x_prime[t]^2));
          y_tilde_real[j,t] <- y_tilde[j,t];
        }
      }
      for (t in 1:T_prime){
        M_tilde[t] <- mean(col(y_tilde_real, t));
        for (j in 1:J_tilde)
          y_tilde_centered[j,t] <- y_tilde_real[j,t] - M_tilde[t];
      }
      Sigma_tilde <- y_tilde_centered'*y_tilde_centered/(J_tilde-1);
      Sigma_tilde <- .5*(Sigma_tilde + Sigma_tilde');
    }
    log_weight <- multi_normal_log(y_prime_bar, M_tilde, Sigma_tilde/J_prime) +
                  multi_normal_log(delta, mu_delta_p, Sigma_delta_p) -
                  multi_normal_log(delta, mu_delta_g, Sigma_delta_g);
    for (k in 1:K)
      output[k] <- delta[k];
    output[K+1] <- log_weight;
    return output;
  }
}
model {
}
