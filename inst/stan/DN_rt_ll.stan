// ll
  vector[Nrate] xbRate  = Xrate * betaRate + offsetInt; 

  for (t in 1:Trate)
    target  += isDependent[t] * xbRate[choseRate[t]] - timespan[t] *
          exp(log_sum_exp(xbRate[startRate[t]:endRate[t]]));

// gq
generated quantities {
  vector[Trate] logLik; 
  {
    vector[Nrate] xbRate = Xrate * betaRate + offsetInt; 

    for (t in 1:Trate)
      logLik[t] = isDependent[t] * xbRate[choseRate[t]] - timespan[t] *
          exp(log_sum_exp(xbRate[startRate[t]:endRate[t]]));
  }
}
