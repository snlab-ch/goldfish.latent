//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  vector[N] y;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real mu;
  real<lower=0> sigma;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  y ~ normal(mu, sigma);
}
//   array[Trate] int<lower=1, upper=kR> zstar;
//   real logp_zstar;
//
//   array[kR] vector[Nrate] xb;
//   for (n in 1:kR)
//     xb[n] = Xrate * beta[n] + offsetInt;
//
//   { // Viterbi algorithm
//     array[Trate, kR] int bpointer; // backpointer to the most likely previous state on the most probable path
//     array[Trate, kR] real delta; // max prob for the sequence up to t
//     // that ends with an emission from state k
//     for(k in 1:kR) // first observation
//       delta[1, k] = log(pi1[k]) + isDependent[1] * xb[n][choseRate[1]] -
//         timespan[1] * exp(log_sum_exp(xb[n][startRate[1]:endRate[1]]));
//
//     for (t in 2:Trate) {
//       for (j in 1:kR) { // i = current (t)
//         delta[t, j] = negative_infinity();
//         for (i in 1:kR) { // i = previous (t-1)
//           real logp;
//           logp = delta[t-1, i] + log(theta[i, j]) +
//             isDependent[t] * xb[n][choseRate[t]] -
//             timespan[t] * exp(log_sum_exp(xb[n][startRate[t]:endRate[t]]));
//             if (logp > delta[t, j]) {
//               bpointer[t, j] = i;
//               delta[t, j] = logp;
//             }
//         }
//       }
//     }
//     logp_zstar = max(delta[Trate]);
//     for (j in 1:kR)
//       if (delta[Trate, j] == logp_zstar)
//         zstar[Trate] = j;
//     for (t in 1:(Trate - 1)) {
//       zstar[Trate - t] = bpointer[Trate - t + 1, zstar[Trate - t + 1]];
//     }
//   }
