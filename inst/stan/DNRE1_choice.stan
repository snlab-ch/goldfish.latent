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
  int Nchoice; // number of events * present actors
  int Tchoice; // number of events
  int Pchoice; // number of covariates

  // covariates, decisions * choice set,
  //   fixed effects (include random effect fix part)
  array[Tchoice] int<lower = 1, upper = Nchoice> choseChoice; // position receiver
  matrix[Nchoice, Pchoice] Xchoice;

  // the starting and ending index observation for each event
  array[Tchoice] int<lower = 1, upper = Nchoice> startChoice;
  array[Tchoice] int<lower = 1, upper = Nchoice> endChoice;

  int A; // number of actors/groups

  int Qchoice; // number of random effects
  matrix[Nchoice, Qchoice] Zchoice;
  // vector[Nchoice] Zchoice;

  array[Nchoice] int<lower = 1, upper = A> senderChoice;


}
parameters {
  vector[Pchoice] betaChoice; // fixed effects, includes average random effect
  // individual random effects
  real<lower=0> sigma; // variance random-effect
  vector[A] gamma_raw; // actor uncentered random-effect
}
// transformed parameters {
  // // sigma in original bugs
  // real<lower=0> sigmasq = square(sigma); // standard deviation of random effect
// }
model {
  // normal prior beta 
  target += std_normal_lpdf(betaChoice);
  
  // priors variance random effects
  target += exponential_lpdf(sigma | 1);
  
  // prior random effects
  target += std_normal_lpdf(gamma_raw);
  
  // ll
  
  // create a temporary holding vector
  vector[Nchoice] xbChoice; 
  {
    vector[A] gamma;
    gamma = sigma * gamma_raw;

    xbChoice = Xchoice * betaChoice +
      to_vector(Zchoice) .* gamma[senderChoice]; 
  }

  for (t in 1:Tchoice)
    target  += xbChoice[choseChoice[t]] -
      log_sum_exp(xbChoice[startChoice[t]:endChoice[t]]);
}
// generated quantities {
//   vector[Tchoice] logLik; 
//   {
//     vector[Nchoice] xbChoice; 
//     vector[A] gamma;
//     gamma = sigma * gamma_raw;
// 
//     xbChoice = Xchoice * betaChoice +
//       to_vector(Zchoice) .* gamma[senderChoice]; 
// 
//     for (t in 1:Tchoice)
//       logLik[t] = xbChoice[choseChoice[t]] -
//         log_sum_exp(xbChoice[startChoice[t]:endChoice[t]]);
//   }
// }
