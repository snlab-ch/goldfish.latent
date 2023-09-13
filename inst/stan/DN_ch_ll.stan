// ll
  vector[Nchoice] xbChoice = Xchoice * betaChoice; 

  for (t in 1:Tchoice)
    target  += xbChoice[choseChoice[t]] -
      log_sum_exp(xbChoice[startChoice[t]:endChoice[t]]);

// gq
generated quantities {
  vector[Tchoice] logLik; 
  {
    vector[Nchoice] xbChoice = Xchoice * betaChoice; 

    for (t in 1:Tchoice)
      logLik[t] = xbChoice[choseChoice[t]] -
        log_sum_exp(xbChoice[startChoice[t]:endChoice[t]]);
  }
}
