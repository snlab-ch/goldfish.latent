//
// DNRE1_choice.stan
// Multinomial choice model with P fixed effect covariates and one random effect
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
  int N_choice; // number of events * present actors
  int T_choice; // number of events
  int P_choice; // number of covariates

  // covariates, decisions * choice set,
  //   fixed effects (include random effect fix part)
  array[T_choice] int<lower = 1, upper = N_choice> chose_choice; // pos receiver
  matrix[N_choice, P_choice] X_choice;

  // the starting and ending index observation for each event
  array[T_choice] int<lower = 1, upper = N_choice> start_choice;
  array[T_choice] int<lower = 1, upper = N_choice> end_choice;

  int A; // number of actors/groups

  int Q_choice; // number of random effects
  matrix[N_choice, Q_choice] Z_choice;
  // vector[N_choice] Z_choice;

  array[N_choice] int<lower = 1, upper = A> sender;
}
transformed data {
  vector[N_choice] Z_temp;
  Z_temp = to_vector(Z_choice);
}
parameters {
  vector[P_choice] beta_choice; // fixed effects, includes average random effect
  // individual random effects
  real<lower=0> sigma; // variance random-effect
  vector[A] gamma_raw; // actor uncentered random-effect
}
transformed parameters {
  vector[A] gamma = sigma * gamma_raw; // Centered random effects
}
model {
  // normal prior beta 
  target += std_normal_lpdf(beta_choice);
  
  // priors variance random effects
  target += exponential_lpdf(sigma | 1);
  
  // prior random effects
  target += std_normal_lpdf(gamma_raw);
  
  // ll
  
  // create a temporary holding vector
  vector[N_choice] xb_choice =
    X_choice * beta_choice + Z_temp .* gamma[sender]; 

  for (t in 1:T_choice)
    target  += xb_choice[chose_choice[t]] -
      log_sum_exp(xb_choice[start_choice[t]:end_choice[t]]);
}
// generated quantities {
//   vector[T_choice] logLik; 
//   {
//     vector[N_choice] xb_choice; 
//     vector[A] gamma;
//     gamma = sigma * gamma_raw;
// 
//     xb_choice = X_choice * beta_choice +
//       to_vector(Z_choice) .* gamma[sender]; 
// 
//     for (t in 1:T_choice)
//       logLik[t] = xb_choice[chose_choice[t]] -
//         log_sum_exp(xb_choice[start_choice[t]:end_choice[t]]);
//   }
// }
