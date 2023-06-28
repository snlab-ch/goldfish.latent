  int Qchoice; // number of random effects
  matrix[Nchoice, Qchoice] Zchoice;

  array[Nchoice] int<lower = 1, upper = A> senderChoice;
