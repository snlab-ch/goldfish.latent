// pm
  real<lower=0> sigma; // variance random-effect
  vector[A] gamma_raw; // actor uncentered random-effect

// tp
transformed parameters {
  // sigma in original bugs
  real<lower=0> sigmasq = square(sigmasq); // standard deviation of random effect
}

// pr:inv-gamma
  target += inv_gamma_lpdf(sigmasq | 1, 1);
  target += std_normal_lpdf(gamma_raw);

// pr:exponential
  target += exponential_lpdf(sigma | 1);
  target += std_normal_lpdf(gamma_raw);
