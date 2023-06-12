// DyNAM_both.stan
// The rate and choice sub-models with intercept model with P covariates
// covariates are standarized and
// the log crude rate is used as offset for the intercept
//
data {
  //
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
  //
  int Nchoice; // number of events * present actors
  int Tchoice; // number of events
  int Pchoice; // number of covariates

  array[Tchoice] int<lower = 1, upper = Nrate> choseChoice; // position sender
  matrix[Nchoice, Pchoice] Xchoice;

  // the starting and ending index observation for each event
  array[Tchoice] int<lower = 1, upper = Nchoice> startChoice;
  array[Tchoice] int<lower = 1, upper = Nchoice> endChoice;
}
parameters {
  vector[Prate] betaRate;
  vector[Pchoice] betaChoice;
}
model {
  // priors
  target += normal_lpdf(betaRate[1] | 0, 5);
  target += std_normal_lpdf(betaRate[2 : ]);
  target += std_normal_lpdf(betaChoice);

  // helper for likelihood
  vector[Trate] logLik;
  {
    vector[Nrate] xbRate;
    xbRate = Xrate * betaRate + offsetInt;

    vector[Nchoice] xbChoice;
    xbChoice = Xchoice * betaChoice;

    int tChoice = 1;
    for (t in 1:Trate) {
      logLik[t] = -timespan[t] *
        exp(log_sum_exp(xbRate[startRate[t]:endRate[t]]));

      if (isDependent[t]) {
        logLik[t] += xbRate[choseRate[t]] +
          xbChoice[choseChoice[tChoice]] -
          log_sum_exp(xbChoice[startChoice[tChoice]:endChoice[tChoice]]);

        tChoice += 1;
      }
    }
  }
  target += logLik;
}
