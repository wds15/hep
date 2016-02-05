data {
  int J;                         // #people in y
  int T;                         // #time points in y
  int K;                         // #parameters in delta
  int K_phi;                     // #parameters in phi
  int y[J,T];
  int J_prime;                   // #people in second dataset
  int T_prime;                   // #time points in second dataset
  int y_prime[J_prime,T_prime];
  vector[T] x;
  vector[T_prime] x_prime;
  vector[K_phi] mu_phi_p;        // prior mean
  cov_matrix[K_phi] Sigma_phi_p; // prior variance
  vector[K] mu_delta_p;          // prior mean
  cov_matrix[K] Sigma_delta_p;   // prior variance
  int n_trials;
  int n_trials_prime;
}
parameters {
  real eta_a[J,K];
  vector[K] mu_a;
  vector<lower=0>[K] sigma_a;    // Now assuming indep
  real b;                        // Shared parameter
  vector[K] delta;
  real eta_prime_a[J_prime,K];
}
transformed parameters {
  real a[J,K];
  vector[K_phi] phi;
  real a_prime[J_prime,K];
  for (j in 1:J)
    for (k in 1:K)
      a[j,k] <- mu_a[k] + sigma_a[k]*eta_a[j,k];
  phi[1] <- mu_a[1];
  phi[2] <- mu_a[2];
  phi[3] <- b;
  phi[4] <- log(sigma_a[1]);
  phi[5] <- log(sigma_a[2]);
  for (j in 1:J_prime)
    for (k in 1:K)
      a_prime[j,k] <- mu_a[k] + delta[k] + sigma_a[k]*eta_prime_a[j,k];
}
model {
  real y_pred[J,T];
  real y_prime_pred[J_prime,T_prime];
  for (j in 1:J)
    for (t in 1:T)
      y_pred[j,t] <- a[j,1] + a[j,2]*x[t] + b*x[t]^2;
  to_array_1d(y) ~ binomial_logit(n_trials, to_array_1d(y_pred));
  to_array_1d(eta_a) ~ normal(0,1);
  phi ~ multi_normal(mu_phi_p, Sigma_phi_p);
  for (j in 1:J_prime)
    for (t in 1:T_prime)
      y_prime_pred[j,t] <- a_prime[j,1] + a_prime[j,2]*x_prime[t] + b*x_prime[t]^2;
  to_array_1d(y_prime) ~ binomial_logit(n_trials_prime, to_array_1d(y_prime_pred));
  to_array_1d(eta_prime_a) ~ normal(0,1);
  delta ~ multi_normal(mu_delta_p, Sigma_delta_p);
}
