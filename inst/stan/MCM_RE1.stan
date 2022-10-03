//
// MCM_RE1.stan
// Multinomial choice model with P fixed effect covariates and one random effect
//
// Update to avoid warnings from deprecated syntax
// non centered parametrization
//
// Subset the vector of probs instead of dot product
// Change to a full mixed formulation,
//   it's easier to decompose effects and add interactions
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

data {
  int N; // number of events * choice set (actors - 1)
  int T; // number of events
  int A; // number of actors
  int P; // number of covariates

  // covariates, decisions * choice set,
  //   fixed effects (include random effect fix part)
  matrix[N, P] X;
  vector[N] Z; // covariate, decisions * choice set, random effect

  array[T] int<lower = 1, upper = N> chose; // position receiver chose

  array[N] int<lower = 1, upper = A> sender; // index for actors

  // the starting and ending index observation for each event
  array[T] int<lower = 1, upper = N> start;
  array[T] int<lower = 1, upper = N> end;
}
parameters {
  vector[P] beta; // fixed effects, includes fixed part random effect
  //real<lower = 0> sigma;
  real<lower=0> sigmasq; // variance of random effect
  vector[A] gamma_raw; // individual uncentered random effects
}
transformed parameters {
  // sigma in original bugs
  real<lower=0> sigma;
  sigma = sqrt(sigmasq); // standard deviation of random effect
}
model {
  // create a temporary holding vector
  vector[N] xb;

  // priors on the parameters
  target += inv_gamma_lpdf(sigmasq | 1, 1);
  //sigma ~ cauchy(0, 5);
  target += normal_lpdf(beta | 0, 4);
  //gamma ~ normal(alpha, sigma);
  target += std_normal_lpdf(gamma_raw);

  // log probabilities of each choice in the dataset
  {
    vector[A] gamma;
    gamma = sigma * gamma_raw;

    xb = X * beta + Z .* gamma[sender];
  }

  for(t in 1:T)
    target  += xb[chose[t]] - log_sum_exp(xb[start[t]:end[t]]);
}
