// Tetraploid F1 test
// Double reduction and no preferential pairing
// parental genotypes are known

functions {
  // alpha Double reduction rate
  // g parent genotype
  // khalf ploidy / 2 + 1
  // return: gamete frequencies of a parent
  vector segfreq4(real alpha, int g) {
    vector[3] p;
    if (g == 0) {
      p[1] = 1.0;
      p[2] = 0.0;
      p[3] = 0.0;
    } else if (g == 1) {
      p[1] = (2.0 + alpha) / 4.0;
      p[2] = 2.0 * (1.0 - alpha) / 4.0;
      p[3] = alpha / 4.0;
    } else if (g == 2) {
      p[1] = (1.0 + 2.0 * alpha) / 6.0;
      p[2] = 4.0 * (1.0 - alpha) / 6.0;
      p[3] = (1.0 + 2.0 * alpha) / 6.0;
    } else if (g == 3) {
      p[1] = alpha / 4.0;
      p[2] = 2.0 * (1.0 - alpha) / 4.0;
      p[3] = (2.0 + alpha) / 4.0;
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
  int<lower=0> x[5]; // genotype counts
  real<lower=0.0,upper=1.0> drbound; // upper bound of double reduction rate
  int<lower=0,upper=4> g1; // first parent genotype
  int<lower=0,upper=4> g2; // second parent genotype
  real<lower=0.0,upper=1.0> mixprop; // mixing component with uniform
}

parameters {
  real<lower=0.0,upper=1.0> tau; // probability of quadrivalent pairing
  real<lower=0.0,upper=drbound> beta; // probability of double reduction given quadrivalent pairing
}

transformed parameters {
  real<lower=0.0,upper=drbound> alpha = tau * beta;
}

model {
  vector[3] p1;
  vector[3] p2;
  vector[5] q;
  vector[5] u = [0.2, 0.2, 0.2, 0.2, 0.2]';
  p1 = segfreq4(alpha, g1);
  p2 = segfreq4(alpha, g2);
  q = [p1[1] * p2[1], p1[1] * p2[2] + p1[2] * p2[1], p1[1] * p2[3] + p1[2] * p2[2] + p1[3] * p2[1], p1[2] * p2[3] + p1[3] * p2[2], p1[3] * p2[3]]';
  q = (1.0 - mixprop) * q + mixprop * u; // mixing to avoid gradient issues
  target += uniform_lpdf(tau | 0.0, 1);
  target += uniform_lpdf(beta | 0.0, drbound);
  target += multinomial_lpmf(x | q);
}
