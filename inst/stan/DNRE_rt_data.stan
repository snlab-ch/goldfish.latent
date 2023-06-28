  int Qrate; // number of random effects
  matrix[Nrate, Qrate] Zrate;

  array[Nrate] int<lower = 1, upper = A> senderRate;
