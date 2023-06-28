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
