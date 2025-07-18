//
// DNRE_Q0_rate.stan
// Rate model with P fixed effect covariates and random intercepts
//
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
  int N_rate; // number of events * present actors
  int T_rate; // number of events
  int P_rate; // number of covariates

  // covariates, decisions * choice set,
  //   fixed effects (include random effect fix part)
  matrix[N_rate, P_rate] X_rate;

  // the starting and ending index observation for each event
  array[T_rate] int<lower = 1, upper = N_rate> start_rate;
  array[T_rate] int<lower = 1, upper = N_rate> end_rate;
  array[T_rate] int<lower = 1, upper = N_rate> chose_rate; // pos receiver

  int A; // number of actors/groups
  array[N_rate] int<lower = 1, upper = A> sender;

  vector<lower=0>[T_rate] timespan;
  array[T_rate] int<lower=0, upper=1> is_dependent;
}
transformed data {
  real log_crude_rate = log(T_rate / (N_rate * mean(timespan)));
}
parameters {
  vector[P_rate] beta_rate; // fixed effects, includes average random effect
  // individual random effects
  real<lower=0> sigma; // variance random-effect
  vector[A] gamma_raw; // actor uncentered random-effect
}
transformed parameters {
  vector[A] gamma = sigma * gamma_raw; // Centered random effects
}
model {
  // normal prior beta
  target += normal_lpdf(beta_rate[1] | log_crude_rate, 4);
  target += std_normal_lpdf(beta_rate[2:]); 
  
  // priors variance random effects
  target += exponential_lpdf(sigma | 1);
  
  // prior random effects
  target += std_normal_lpdf(gamma_raw);
  
  // ll
  
  // create a temporary holding vector
  vector[N_rate] xb_rate =
    X_rate * beta_rate + gamma[sender]; 

  for (t in 1:T_rate)
    if (timespan[t] > 0)
      target  += is_dependent[t] * xb_rate[chose_rate[t]] -
        timespan[t] * exp(log_sum_exp(xb_rate[start_rate[t]:end_rate[t]]));
}
