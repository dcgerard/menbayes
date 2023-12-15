// Tetraploid F1 test
// Double reduction and preferential pairing
// neither offspring nor parental genotypes are known.

functions {
  // alpha Double reduction rate
  // xi Preferential pairing rate
  // g parent genotype
  // khalf ploidy / 2 + 1
  // return: gamete frequencies of a parent
  vector segfreq4(real alpha, real xi, int g) {
    vector[3] p;
    if (g == 0) {
      p[1] = 1.0;
      p[2] = 0.0;
      p[3] = 0.0;
    } else if (g == 1) {
      p[1] = 0.5 + 0.25 * alpha;
      p[2] = 0.5 - 0.5 * alpha;
      p[3] = 0.25 * alpha;
    } else if (g == 2) {
      p[1] = 0.5 * alpha + 0.25 * (1.0 - alpha) * (1.0 - xi);
      p[2] = 0.5 * (1.0 - alpha) * (1.0 + xi);
      p[3] = 0.5 * alpha + 0.25 * (1.0 - alpha) * (1.0 - xi);
    } else if (g == 3) {
      p[1] = 0.25 * alpha;
      p[2] = 0.5 - 0.5 * alpha;
      p[3] = 0.5 + 0.25 * alpha;
    } else if (g == 4) {
      p[1] = 0.0;
      p[2] = 0.0;
      p[3] = 1.0;
    } else {
      // do nothing
    }
    return p;
  }
}

data {
  int N;
  vector[5] p1_gl; // genotype log-likelihoods for parent 1
  vector[5] p2_gl; // genotype log-likelihoods for parent 2
  matrix[N, 5] gl; // genotype log-likelihoods for offspring
  real<lower=0.0,upper=1.0> drbound; // upper bound of double reduction rate
  real<lower=0.0,upper=1.0> mixprop; // mixing component with uniform
  real<lower=0.0> shape1; // shape 1 of beta
  real<lower=0.0> shape2; // shape 2 of beta
}

parameters {
  real<lower=0.0,upper=drbound> beta; // double reduction given quad
  real<lower=0.0,upper=1.0> tau; // prob quad
  real<lower=0.0,upper=1.0> gamma1; // prob AA_aa
  real<lower=0.0,upper=1.0> gamma2; // prob AA_aa
}

transformed parameters {
  real<lower=0.0,upper=drbound> alpha = beta * tau; // double reduction rate
  real<lower=0.0,upper=1.0> eta = (1.0 - beta) * tau / ((1.0 - beta) * tau + (1.0 - tau)); // prob quad given no dr
  real<lower=0.0,upper=1.0> xi1 = eta / 3.0 + (1.0 - eta) * gamma1; // preferential pairing rate
  real<lower=0.0,upper=1.0> xi2 = eta / 3.0 + (1.0 - eta) * gamma2; // preferential pairing rate
  matrix[5, 5] glmat;
  for (i in 1:5) {
    vector[3] p1;
    p1 = segfreq4(alpha, xi1, i - 1);
    for (j in 1:5) {
      vector[5] u = [0.2, 0.2, 0.2, 0.2, 0.2]';
      vector[3] p2;
      vector[5] q;
      p2 = segfreq4(alpha, xi2, j - 1);
      q = [p1[1] * p2[1], p1[1] * p2[2] + p1[2] * p2[1], p1[1] * p2[3] + p1[2] * p2[2] + p1[3] * p2[1], p1[2] * p2[3] + p1[3] * p2[2], p1[3] * p2[3]]';
      q = (1.0 - mixprop) * q + mixprop * u; // mixing to avoid gradient issues
      glmat[i, j] = p1_gl[i] + p2_gl[j];
      for (ind in 1:N) {
        glmat[i, j] += log_sum_exp(to_vector(gl[ind]) + log(q));
      }
    }
  }
}

model {
  target += uniform_lpdf(beta | 0.0, drbound);
  target += uniform_lpdf(tau | 0.0, 1.0);
  target += beta_lpdf(gamma1 | shape1, shape2);
  target += beta_lpdf(gamma2 | shape1, shape2);
  target += log_sum_exp(glmat);
}
