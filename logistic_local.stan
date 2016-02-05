data {
  int J;                         // #people in y
  int T;                         // #time points in y
  int K;                         // #parameters in delta
  int K_phi;                     // #parameters in phi
  int y[J,T];
  vector[T] x;
  vector[K_phi] mu_phi_p;        // prior mean
  cov_matrix[K_phi] Sigma_phi_p; // prior variance
  vector[K_phi] mu_phi_g;        // pseudo-prior mean
  cov_matrix[K_phi] Sigma_phi_g; // pseudo-prior variance
  int n_trials;                  // binomial "n"
}
parameters {
  real eta_a[J,K];
  real<lower=0> sigma_y;
  vector[K] mu_a;
  vector<lower=0>[K] sigma_a;    // Now assuming indep
  real b;                        // Shared parameter
}
transformed parameters {
  real a[J,K];
  vector[K_phi] phi;
  for (j in 1:J)
    for (k in 1:K)
      a[j,k] <- mu_a[k] + sigma_a[k]*eta_a[j,k];
  phi[1] <- mu_a[1];
  phi[2] <- mu_a[2];
  phi[3] <- b;
  phi[4] <- log(sigma_a[1]);
  phi[5] <- log(sigma_a[2]);
}
model {
  real y_pred[J,T];
  for (j in 1:J)
    for (t in 1:T)
      y_pred[j,t] <- a[j,1] + a[j,2]*x[t] + b*x[t]^2;
  to_array_1d(y) ~ binomial_logit(n_trials, to_array_1d(y_pred));
  to_array_1d(eta_a) ~ normal(0,1);
  phi ~ multi_normal(mu_phi_g, Sigma_phi_g);
}
generated quantities {
  real log_weight;
  log_weight <- multi_normal_log(phi, mu_phi_p, Sigma_phi_p) -
                multi_normal_log(phi, mu_phi_g, Sigma_phi_g);
}
