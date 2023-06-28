// parms
  real<lower=0> sigmasq; // variance random-effect
  vector[A] gamma_raw; // actor uncentered random-effect 
// transformedParms
transformed parameters {
  // sigma in original bugs
  real<lower=0> sigma;
  sigma = sqrt(sigmasq); // standard deviation of random effect
}
// prior
  target += inv_gamma_lpdf(sigmasq | 1, 1);
  target += std_normal_lpdf(gamma_raw);
