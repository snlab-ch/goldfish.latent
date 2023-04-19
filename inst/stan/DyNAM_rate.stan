// DyNAM_rate.stan
// The rate model with intercept model with P covariates
//
data {
  int Nrate; // number of events * present actors
  int Trate; // number of events
  int Prate; // number of covariates

  array[Trate] int<lower = 0, upper = Nrate> choseRate; // position sender
  matrix[Nrate, Prate] Xrate;

  // the starting and ending index observation for each event
  array[Trate] int<lower = 1, upper = Nrate> startRate;
  array[Trate] int<lower = 1, upper = Nrate> endRate;

  //
  array[Trate] real<lower = 0> timespan;
  array[Trate] int<lower = 0, upper = 1> isDependent;

  real offsetInt;
}
parameters {
  vector[Prate] beta;
}
model {
  // priors
  target += normal_lpdf(beta| 0, 4);

  // helper for likelihood
  {
    vector[Nrate] xb;
    xb = Xrate * beta + offsetInt;

    for(t in 1:Trate) {
      if (timespan[t] > 0)
        target += isDependent[t] * xb[choseRate[t]] -
          timespan[t] * exp(log_sum_exp(xb[startRate[t]:endRate[t]]));
    }
  }
}
