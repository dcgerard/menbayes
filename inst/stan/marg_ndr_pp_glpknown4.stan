// Tetraploid F1 test
// No double reduction and preferential pairing
// parental genotypes are known, offspring genotype likelihoods are used

functions {
  // xi Preferential pairing rate
  // g parent genotype
  // khalf ploidy / 2 + 1
  // return: gamete frequencies of a parent
  vector segfreq4(real xi, int g) {
    vector[3] p;
    if (g == 0) {
      p[1] = 1.0;
      p[2] = 0.0;
      p[3] = 0.0;
    } else if (g == 1) {
      p[1] = 0.5;
      p[2] = 0.5;
      p[3] = 0.0;
    } else if (g == 2) {
      p[1] = 0.25 * (1.0 - xi);
      p[2] = 0.5 * (1.0 + xi);
      p[3] = 0.25 * (1.0 - xi);
    } else if (g == 3) {
      p[1] = 0.0;
      p[2] = 0.5;
      p[3] = 0.5;
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
  matrix[N, 5] gl; // genotype log-likelihoods for offspring
  int<lower=0,upper=4> g1; // first parent genotype
  int<lower=0,upper=4> g2; // second parent genotype
  real<lower=0.0,upper=1.0> mixprop; // mixing component with uniform
}

parameters {
  real<lower=0,upper=1> xi; // preferential pairing rate
}

model {
  vector[3] p1;
  vector[3] p2;
  vector[5] q;
  vector[5] u = [0.2, 0.2, 0.2, 0.2, 0.2]';
  p1 = segfreq4(xi, g1);
  p2 = segfreq4(xi, g2);
  q = [p1[1] * p2[1], p1[1] * p2[2] + p1[2] * p2[1], p1[1] * p2[3] + p1[2] * p2[2] + p1[3] * p2[1], p1[2] * p2[3] + p1[3] * p2[2], p1[3] * p2[3]]';
  q = (1.0 - mixprop) * q + mixprop * u; // mixing to avoid gradient issues
  target += beta_lpdf(xi | 1.0, 2.0);
  for (ind in 1:N) {
    target += log_sum_exp(to_vector(gl[ind]) + log(q));
  }
}