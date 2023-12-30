// Tetraploid F1 test
// Double reduction and preferential pairing
// parental genotypes are known

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
  int<lower=0> x[5]; // genotype counts
  int<lower=0,upper=4> g1; // first parent genotype
  int<lower=0,upper=4> g2; // second parent genotype
  real<lower=0.0,upper=1.0> mixprop; // mixing component with uniform
  real<lower=0.0> shape1; // shape 1 of beta
  real<lower=0.0> shape2; // shape 2 of beta
  real<lower=0.0,upper=1.0> alpha; // known double reduction rate
}

parameters {
  real<lower=0.0,upper=1.0> gamma1; // preferential pairing rate
  real<lower=0.0,upper=1.0> gamma2; // preferential pairing rate
}

model {
  vector[3] p1;
  vector[3] p2;
  vector[5] q;
  vector[5] u = [0.2, 0.2, 0.2, 0.2, 0.2]';
  p1 = segfreq4(alpha, gamma1, g1);
  p2 = segfreq4(alpha, gamma2, g2);
  q = [p1[1] * p2[1], p1[1] * p2[2] + p1[2] * p2[1], p1[1] * p2[3] + p1[2] * p2[2] + p1[3] * p2[1], p1[2] * p2[3] + p1[3] * p2[2], p1[3] * p2[3]]';
  q = (1.0 - mixprop) * q + mixprop * u; // mixing to avoid gradient issues
  target += beta_lpdf(gamma1 | shape1, shape2);
  target += beta_lpdf(gamma2 | shape1, shape2);
  target += multinomial_lpmf(x | q);
}
