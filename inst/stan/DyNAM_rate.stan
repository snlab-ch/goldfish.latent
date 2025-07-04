// DyNAM_rate.stan
// The rate model with intercept model with P covariates
//
data {
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

  //real offsetInt;
}
transformed data {
  // Trate / sum(timespan) / mean(actors)
  real log_crude_rate = log(T_rate / (N_rate * mean(timespan)));
}
parameters {
  vector[P_rate] beta_rate;
}
model {
  // priors
  target += normal_lpdf(beta_rate[1] | 0, 4);
  target += std_normal_lpdf(beta_rate[2:]);

  // helper for likelihood
  {
    vector[N_rate] xb;
    xb = X_rate * beta_rate + log_crude_rate;

    for(t in 1:T_rate) {
      if (timespan[t] > 0)
        target += (is_dependent[t] ? xb[chose_rate[t]] : 0) -
          timespan[t] * exp(log_sum_exp(xb[start_rate[t]:end_rate[t]]));
    }
  }
}
