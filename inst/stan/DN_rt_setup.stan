// dt
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
// pm
  vector[Prate] betaRate;
// pr:normal
  target += normal_lpdf(betaRate[1] | 0, 10);
  target += std_normal_lpdf(betaRate[2: ]);
// pr:t-student
  target += student_t_lpdf(betaRate[1] | 3, 0, 10);
  target += student_t_lpdf(betaRate[2: ] | 3, 0, 1);
