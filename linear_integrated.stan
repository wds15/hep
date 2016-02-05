data {
  int J;                         // #people in y
  int T;                         // #time points in y
  int K;                         // #parameters in delta
  int K_phi;                     // #parameters in phi
  matrix[J,T] y;
  int J_prime;                   // #people in second dataset
  int T_prime;                   // #time points in second dataset
  vector[T_prime] y_prime_bar;
  vector[T] x;
  vector[T_prime] x_prime;
  vector[K_phi] mu_phi_p;        // prior mean
  cov_matrix[K_phi] Sigma_phi_p; // prior variance
  vector[K] mu_delta_p;          // prior mean
  cov_matrix[K] Sigma_delta_p;   // prior variance
}
transformed data {
  vector[T_prime] ones;
  matrix[T_prime,K] X_prime;
  cov_matrix[T_prime] identity;
  ones <- rep_vector(1, T_prime);
  X_prime <- append_col(ones, x_prime);
  identity <- diag_matrix(ones);
}
parameters {
  matrix[J,K] eta_a;
  real<lower=0> sigma_y;
  vector[K] mu_a;
  vector<lower=0>[K] sigma_a;    // Now assuming indep
  real b;                        // Shared parameter
  vector[K] delta;
  matrix[J_prime,K] eta_prime_a;
}
transformed parameters {
  matrix[J,K] a;
  vector[K_phi] phi;
  for (j in 1:J)
    for (k in 1:K)
      a[j,k] <- mu_a[k] + sigma_a[k]*eta_a[j,k];
  phi[1] <- mu_a[1];
  phi[2] <- mu_a[2];
  phi[3] <- b;
  phi[4] <- log(sigma_a[1]);
  phi[5] <- log(sigma_a[2]);
  phi[6] <- log(sigma_y);
}
model {
  matrix[J,T] y_pred;
  vector[T_prime] y_prime_pred;
  matrix[T_prime,T_prime] Sigma_y_prime_bar;
  for (j in 1:J)
    for (t in 1:T)
      y_pred[j,t] <- a[j,1] + a[j,2]*x[t] + b*x[t]^2;
  to_vector(y) ~ normal(to_vector(y_pred), sigma_y);
  to_vector(eta_a) ~ normal(0,1);
  Sigma_y_prime_bar <- (X_prime*sigma_a*sigma_a'*X_prime' + sigma_y^2*identity)/J_prime;
  for (t in 1:T_prime)
    y_prime_pred[t] <- mu_a[1] + delta[1] + (mu_a[2] + delta[2])*x_prime[t] + b*x_prime[t]^2;
  y_prime_bar ~ multi_normal(y_prime_pred, Sigma_y_prime_bar);
  phi ~ multi_normal(mu_phi_p, Sigma_phi_p);
  delta ~ multi_normal(mu_delta_p, Sigma_delta_p);
}