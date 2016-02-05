functions {

  // takes sampling parameters and produces model prediction for all
  // patients
  matrix evaluate_model(vector x, int[] DRUG,
			real lalpha_0, real lalpha_s, real lkappa, real lEmax_s, real delta, vector eta_lalpha_0, vector eta_lalpha_s, real sigma_lalpha_0, real sigma_lalpha_s) {
    matrix[num_elements(eta_lalpha_0), num_elements(x)] lpred;
    vector[num_elements(eta_lalpha_0)] lalpha_0_j;
    vector[num_elements(eta_lalpha_0)] lalpha_s_j;
    real kout;
    vector[3] FlEmax_s;
    
    lalpha_0_j  <- lalpha_0 + sigma_lalpha_0 * eta_lalpha_0;
    lalpha_s_j  <- lalpha_s + sigma_lalpha_s * eta_lalpha_s;
    kout <- exp((lkappa - lalpha_s)/2.);

    FlEmax_s[1] <- 0.;
    FlEmax_s[2] <- log1p_exp(lEmax_s);    // =  log( 1+exp(lEmax_s) )
    FlEmax_s[3] <- log1p_exp(lEmax_s + delta);
    
    // a + (b-a) * exp(-kout*x) = a * (1 - exp(-kout*x)) + b * exp(-kout*x)
    for(j in 1:num_elements(eta_lalpha_s)) {
      for(t in 1:num_elements(x))
	if(fabs(x[t]) < 1e-8) {
	  // for x == 0, then this log formulation becomes numerically
	  // unstable
	  lpred[j,t] <-  lalpha_0_j[j];
	} else {
	  lpred[j,t] <- log_sum_exp(lalpha_s_j[j] + FlEmax_s[DRUG[j]+1] + log1m_exp(-kout * x[t]),
				    lalpha_0_j[j] - kout * x[t]);
	}
    }

    return(lpred);
  }

  void pretty_print(matrix x) {
    if (rows(x) == 0) {
      print("empty matrix");
      return;
    }
    for (m in 1:rows(x)) {
      row_vector[cols(x)] rv;
      for (n in 1:cols(x))
	rv[n] <- round(1000*x[m,n])/1000.;
      print("row ", m, " = ", rv);      
    }
  }

  // numerically robust covariance estimate
  matrix robust_cov(matrix y) {
    matrix[cols(y),cols(y)] cov;
    matrix[rows(y),cols(y)] yc;
    
    // center column wise
    for (i in 1:cols(y)) {
      real m;
      m <- mean(col(y, i));
      for(j in 1:rows(y)) {
	yc[j,i] <- y[j,i] - m;
      }
    }
    
    cov <- yc' * yc /(rows(y) - 1);
    return(0.5 * (cov + cov'));
  }

  vector colMeans(matrix y) {
    vector[cols(y)] m;
    for (i in 1:cols(y))
      m[i] <- mean(col(y, i));
    return(m);
  }
      

  vector pkpd_impsampling_rng(int J_tilde, int J_prime, int K, vector y_prime_bar, vector x_prime, int DRUG_prime,
			      vector mu_delta_p, matrix Sigma_delta_p, vector mu_delta_g, matrix Sigma_delta_g, 
			      real lalpha_0, real lalpha_s, real lkappa, real lEmax_s, real sigma_lalpha_0, real sigma_lalpha_s, real sigma_y) {
    vector[num_elements(x_prime)] y_prime_mean;
    matrix[num_elements(x_prime),num_elements(x_prime)] y_prime_var;
    vector[J_tilde] eta_lalpha_0_sim;
    vector[J_tilde] eta_lalpha_s_sim;
    int DRUG[J_tilde];
    matrix[J_tilde,num_elements(x_prime)] y_prime_sim;
    vector[K+1] output;
    vector[K] delta;
    real log_weight;
    real log_dens_y_prime_bar;
    real log_dens_delta_p;
    real log_dens_delta_g;
    delta <- multi_normal_rng(mu_delta_g, Sigma_delta_g);
    for(j in 1:J_tilde) {
      eta_lalpha_0_sim[j]  <- normal_rng(0, 1);
      eta_lalpha_s_sim[j]  <- normal_rng(0, 1);
      DRUG[j] <- DRUG_prime;
    }
    y_prime_sim <- evaluate_model(x_prime, DRUG,
				  lalpha_0, lalpha_s, lkappa, lEmax_s, delta[1],
				  eta_lalpha_0_sim, eta_lalpha_s_sim,
				  sigma_lalpha_0, sigma_lalpha_s);

    y_prime_mean <- colMeans(y_prime_sim);
    y_prime_var <- robust_cov(y_prime_sim);
    
    // add in additive noise term into covariance matrix
    y_prime_var <- y_prime_var + diag_matrix(rep_vector(sigma_y^2, num_elements(x_prime)));

    log_dens_y_prime_bar <- multi_normal_log(log(y_prime_bar), y_prime_mean, y_prime_var/J_prime);
    log_dens_delta_p <- multi_normal_log(delta, mu_delta_p, Sigma_delta_p);
    log_dens_delta_g <- multi_normal_log(delta, mu_delta_g, Sigma_delta_g);
    log_weight <- log_dens_y_prime_bar + log_dens_delta_p - log_dens_delta_g;
    for (k in 1:K)
      output[k] <- delta[k];
    output[K+1] <- log_weight;

    return output;
  }

}
data {
  int<lower=1> K;                         // #parameters in delta
  int<lower=K> K_phi;                     // #parameters in phi
  int<lower=0> J;
  int<lower=0> T;
  int<lower=0> J_prime;
  int<lower=0> T_prime;
  matrix[J,T] y;
  matrix[J_prime,T_prime] y_prime;
  vector[T] x;
  vector[T_prime] x_prime;
  int<lower=0> J_tilde;
  int<lower=0,upper=2> DRUG[J+J_prime];      // drug administered (0 for Placebo, 1 for directly fit and 2 for fitting via summaries)
  int<lower=0,upper=1> fit_all;
  
  vector[K]     mu_delta_p;        // prior delta mean
  vector[K_phi] mu_phi_p;          // prior mean
  vector[K_phi] mu_phi_g;          // pseudo-prior mean
  cov_matrix[K]     Sigma_delta_p; // prior delta variance
  cov_matrix[K_phi] Sigma_phi_p;   // prior variance
  cov_matrix[K_phi] Sigma_phi_g;   // pseudo-prior variance
}
transformed data {
  vector[T_prime] y_prime_bar;
  matrix[J,T] ly;
  matrix[J_prime,T_prime] ly_prime;
  int J_tot;
  cholesky_factor_cov[K]     L_Sigma_delta_p;
  cholesky_factor_cov[K_phi] L_Sigma_phi_p;
  cholesky_factor_cov[K_phi] L_Sigma_phi_g;

  L_Sigma_delta_p <- cholesky_decompose(Sigma_delta_p);
  L_Sigma_phi_p   <- cholesky_decompose(Sigma_phi_p);
  L_Sigma_phi_g   <- cholesky_decompose(Sigma_phi_g);

  J_tot <- J + J_prime;

  ly <- log(y);
  ly_prime <- log(y_prime);

  // we use geometric means
  y_prime_bar <- exp(colMeans(ly_prime));
}
parameters {
  // population parameters
  real lalpha_0;              // log-BCVA latent at baseline
  real lalpha_s;              // alpha_s = kin_sp/kout
  real lkappa;                // kout = 1/exp(ltau)
  real lEmax_s;               // emax in ss for reference drug
  vector[K] delta;            // change in Emax for second drug

  real<lower=0> sigma_lalpha_0;
  real<lower=0> sigma_lalpha_s;

  vector[J_tot] eta_lalpha_0;
  vector[J_tot] eta_lalpha_s;

  real<lower=0> sigma_y;            // prediction noise term
}
transformed parameters {
  real y_log;
  real y_prime_log;
  matrix[J,T] y_pred;
  matrix[J_prime,T_prime] y_prime_pred;
  vector[K_phi] phi;

  y_pred       <- evaluate_model(x, segment(DRUG, 1, J),
				 lalpha_0, lalpha_s, lkappa, lEmax_s, delta[1],
				 segment(eta_lalpha_0, 1, J), segment(eta_lalpha_s, 1, J),
				 sigma_lalpha_0, sigma_lalpha_s);

  y_prime_pred <- evaluate_model(x_prime, segment(DRUG, J+1, J_prime),
				 lalpha_0, lalpha_s, lkappa, lEmax_s, delta[1],
				 segment(eta_lalpha_0, J+1, J_prime), segment(eta_lalpha_s, J+1, J_prime), 
				 sigma_lalpha_0, sigma_lalpha_s);

  y_log       <- normal_log(to_vector(ly),       to_vector(y_pred),       sigma_y);
  y_prime_log <- normal_log(to_vector(ly_prime), to_vector(y_prime_pred), sigma_y);

  phi[1] <- lalpha_0;
  phi[2] <- lalpha_s;
  phi[3] <- lkappa;
  phi[4] <- lEmax_s;
  phi[5] <- log(sigma_lalpha_0);
  phi[6] <- log(sigma_lalpha_s);
  phi[7] <- log(sigma_y);
}
model {
  increment_log_prob(y_log);

  // we only add the data from the second data set if requested
  if(fit_all == 1)
    increment_log_prob(y_prime_log);

  // prior for the delta of the second drug
  delta ~ multi_normal_cholesky(mu_delta_p, L_Sigma_delta_p);

  eta_lalpha_0  ~ normal(0, 1);
  eta_lalpha_s  ~ normal(0, 1);

  phi ~ multi_normal_cholesky(mu_phi_g, L_Sigma_phi_g);
}
generated quantities {
  real log_weight;
  log_weight <- multi_normal_cholesky_log(phi, mu_phi_p, L_Sigma_phi_p) -
                multi_normal_cholesky_log(phi, mu_phi_g, L_Sigma_phi_g);
}