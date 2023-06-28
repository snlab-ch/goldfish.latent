// loglikelihood
  vector[Nchoice] xbChoice; 
  {
    vector[A] gamma;
    gamma = sigma * gamma_raw;

    xbChoice = Xchoice * betaChoice +
      to_vector(col(Zchoice, 1)) .* gamma[senderChoice]; 
  }

  for(t in 1:Tchoice)
    target  += xbChoice[choseChoice[t]] -
      log_sum_exp(xbChoice[startChoice[t]:endChoice[t]]);
// generateQ
generated quantities {
  vector[Tchoice] logLik; 
  {
    vector[Nchoice] xbChoice; 
    vector[A] gamma;
    gamma = sigma * gamma_raw;

    xbChoice = Xchoice * betaChoice +
      to_vector(col(Zchoice, 1)) .* gamma[senderChoice]; 

    for(t in 1:Tchoice)
      logLik[t] = xbChoice[choseChoice[t]] -
        log_sum_exp(xbChoice[startChoice[t]:endChoice[t]]);
  }
}
