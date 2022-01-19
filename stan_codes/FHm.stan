data {
  int<lower=0>            m2;   // Number of sampled small areas
  int<lower=0>            m1;   // Number of non-sampled small areas
  int<lower=0>             p;   // Number of auxiliary variables
  real                 y[m2];   // Direct Estimate
  real<lower=0>      sDi[m2];   // Sampling error
  vector[p]         X[m1+m2];   // Auxiliary variable
  matrix[m2,m1+m2]          M;  // M*theta = theta_(2)
}
parameters {
  real<lower=0>  sigma_sq;    // Model error variance
  vector[p]          beta;      // Model regression parameter
  vector[m1+m2]      theta;  // Area characteristics
}
transformed parameters {
  vector[m2]       Mtheta;
  real             mu[m1+m2];
  for( i in 1:(m1+m2) ){
    mu[i] = X[i]'*beta;
    }
  Mtheta = M*theta;
}
model {
  for( i in 1:(m1+m2)){
  theta[i] ~ normal( mu[i], sqrt(sigma_sq) );
  }
  for( i in 1:m2){
  y[i] ~ normal( Mtheta[i], sDi[i]);
  }
}
generated quantities {
  vector[m2] log_lik;
  for (n in 1:m2) {
    log_lik[n] = normal_lpdf(y[n] | Mtheta[n], sDi[n]);
  }
}
