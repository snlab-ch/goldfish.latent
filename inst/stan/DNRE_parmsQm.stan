// parms
  cholesky_factor_corr[Q] L_Omega; // correlation RE
  vector<lower=0>[Q] L_sigma; // scale RE  real<lower=0> 

  array[A] vector[Q] gamma_raw;  // group Random Effect
// transformedParms
transformed parameters {
  // Sigma matrix
  matrix[2, 2] Sigma = diag_pre_multiply(L_sigma, L_Omega) *
    diag_pre_multiply(L_sigma, L_Omega)';
}
// prior
  target += lkj_corr_cholesky_lpdf(L_Omega | 2); // bigger than one
  target += exponential_lpdf(L_sigma | 1); // prior over sd, instead of cauchy
