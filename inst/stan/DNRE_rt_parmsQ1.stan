// loglikelihood
  vector[Nrate] xbRate; 
  {
    vector[A] gamma;
    gamma = sigma * gamma_raw;

    xbRate = Xrate * betaRate + offsetInt +
      to_vector(col(Zrate, 1)) .* gamma[senderRate]; 
  }

  for(t in 1:Trate)
    target  += isDependent[t] * xbRate[choseRate[t]] - timespan[t] *
          exp(log_sum_exp(xbRate[startRate[t]:endRate[t]]));

// generateQ
generated quantities {
  vector[Trate] logLik; 
  {
    vector[Nrate] xbRate; 
    vector[A] gamma;
    gamma = sigma * gamma_raw;

    xbRate = Xrate * betaRate + offsetInt +
      to_vector(col(Zrate, 1)) .* gamma[senderRate]; 

    for(t in 1:Trate)
      logLik[t] = isDependent[t] * xbRate[choseRate[t]] - timespan[t] *
          exp(log_sum_exp(xbRate[startRate[t]:endRate[t]]));
  }
}
