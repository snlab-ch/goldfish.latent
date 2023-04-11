// DyNAM_choice_poisson.stan
// Multinomial choice model with P covariates represented as a poisson glm
//
data {
  int N; // number of events * choice set (actors - 1)
  int T; // number of events
  int P; // number of covariates

  array[N] int<lower = 0> selected; // dummy receiver selected
  array[N] int<lower = 0, upper = T> event; // event id
  matrix[N, P] X;
}
parameters {
  vector[P] beta;
  //real<lower=0> sigma_alpha; // variance of random intercepts
  vector[T] alpha_raw;
}
model {
  // priors
  target += normal_lpdf(beta | 0, 4);
  //target += exponential_lpdf(sigma_alpha | 1);
  target += std_normal_lpdf(alpha_raw);

  // helper for likelihood
  vector[N] alpha = alpha_raw[event] * 10;

  target += poisson_log_glm_lpmf(selected | X, alpha, beta);
}
