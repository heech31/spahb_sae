data {
  int<lower=0>       m;    // Number of small areas
  int<lower=0>       p;    // Number of auxiliary variables
  real            y[m];    // Direct Estimate
  real<lower=0> sDi[m];    // Sampling error
  vector[p]        X[m];    // Auxiliary variable
}
parameters {
  real<lower=0>    sigma_sq;     // Log model error variance
  vector[p]          beta;     // Model regression parameter
  real           theta[m];     // Area characteristics
}
transformed parameters {
  real              mu[m];
  for( i in 1:m){
    mu[i] = X[i]'*beta;
    }
}
model {
  for( i in 1:m){
  theta[i] ~ normal( mu[i], sqrt(sigma_sq) );
  }
  for( i in 1:m){
  y[i] ~ normal( theta[i], sDi[i]);
  }
}
