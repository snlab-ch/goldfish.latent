// DyNAM_both.stan
// The rate and choice sub-models with intercept model with P covariates
// covariates are standarized and
// the log crude rate is used as offset for the intercept
//
data {
  //
  int N_rate; // number of events * present actors
  int T_rate; // number of events
  int P_rate; // number of covariates

  array[T_rate] int<lower = 0, upper = N_rate> chose_rate; // position sender
  matrix[N_rate, P_rate] X_rate;

  // the starting and ending index observation for each event
  array[T_rate] int<lower = 1, upper = N_rate> start_rate;
  array[T_rate] int<lower = 1, upper = N_rate> end_rate;

  //
  array[T_rate] real<lower = 0> timespan;
  array[T_rate] int<lower = 0, upper = 1> is_dependent;

  real offset_int;
  //
  int N_choice; // number of events * present actors
  int T_choice; // number of events
  int P_choice; // number of covariates

  array[T_choice] int<lower = 1, upper = N_choice> chose_choice; // pos receiver
  matrix[N_choice, P_choice] X_choice;

  // the starting and ending index observation for each event
  array[T_choice] int<lower = 1, upper = N_choice> start_choice;
  array[T_choice] int<lower = 1, upper = N_choice> end_choice;
}
parameters {
  vector[P_rate] beta_rate;
  vector[P_choice] beta_choice;
}
model {
  // priors
  target += normal_lpdf(beta_rate[1] | 0, 5);
  target += std_normal_lpdf(beta_rate[2 : ]);
  target += std_normal_lpdf(beta_choice);

  // helper for likelihood
  vector[T_rate] logLik;
  {
    vector[N_rate] xb_rate;
    xb_rate = X_rate * beta_rate + offset_int;

    vector[N_choice] xb_choice;
    xb_choice = X_choice * beta_choice;

    int t_choice = 1;
    for (t in 1:T_rate) {
      logLik[t] = -timespan[t] *
        exp(log_sum_exp(xb_rate[start_rate[t]:end_rate[t]]));

      if (is_dependent[t]) {
        logLik[t] += xb_rate[chose_rate[t]] +
          xb_choice[chose_choice[t_choice]] -
          log_sum_exp(xb_choice[start_choice[t_choice]:end_choice[t_choice]]);

        t_choice += 1;
      }
    }
  }
  target += logLik;
}
