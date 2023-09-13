// loglikelihood
  vector[Nrate] xbRate; 
  {
    vector[A] gamma;
    gamma = sigma * gamma_raw;

    xbRate = Xrate * betaRate +
      to_vector(col(Zrate, 1)) .* gamma[senderRate]; 
  }

  for(t in 1:Trate)
    target  += xbRate[choseRate[t]] -
      log_sum_exp(xbRate[startRate[t]:endRate[t]]);
// generateQ
generated quantities {
  vector[Trate] logLik; 
  {
    vector[Nrate] xbRate; 
    vector[A] gamma;
    gamma = sigma * gamma_raw;

    xbRate = Xrate * betaRate +
      to_vector(col(Zrate, 1)) .* gamma[senderRate]; 

    for(t in 1:Trate)
      logLik[t] = xbRate[choseRate[t]] -
        log_sum_exp(xbRate[startRate[t]:endRate[t]]);
  }
}
