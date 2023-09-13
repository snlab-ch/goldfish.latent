// dt
  int Nchoice; // number of events * present actors
  int Tchoice; // number of events
  int Pchoice; // number of covariates

  array[Tchoice] int<lower = 1, upper = Nchoice> choseChoice; // position receiver
  matrix[Nchoice, Pchoice] Xchoice;

  // the starting and ending index observation for each event
  array[Tchoice] int<lower = 1, upper = Nchoice> startChoice;
  array[Tchoice] int<lower = 1, upper = Nchoice> endChoice;
// pm
  vector[Pchoice] betaChoice;
// pr:normal
  target += std_normal_lpdf(betaChoice);
// pr:t-student
  target += student_t_lpdf(betaChoice | 3, 0, 1);