// Polyploid F1 test

functions {
  // p1 gamete frequencies from parent 1.
  // p2 gamete frequencies from parent 2.
  // K ploidy
  // khalf K / 2 + 1 so stan does not complain about integer division
  vector convolve(vector p1, vector p2, int K, int khalf) {
    vector[K+1] q;
    for (k in 1:(K+1)) {
      int iup = min(k - 1, khalf - 1);
      int ilo = max(0, k - khalf);
      q[k] = 0.0;
      for (i in ilo:iup) {
        q[k] += p1[i + 1] * p2[k - i];
      }
    }
    return q;
  }
}

data {
  int<lower=2> ploidy;
  int<lower=1> phalf; // ploidy / 2 + 1
  int<lower=0> x[ploidy + 1]; // genotype counts
  real<lower=0,upper=1> mixprop;
  vector[phalf] beta; // dirichlet concentration prior.
}

parameters {
  simplex[phalf] p1;
  simplex[phalf] p2;
}

model {
  vector[ploidy + 1] q;
  vector[ploidy + 1] u = rep_vector(1 / (ploidy * 1.0 + 1), ploidy + 1);
  q = convolve(p1, p2, ploidy, phalf);
  q = (1.0 - mixprop) * q + mixprop * u;
  target += dirichlet_lpdf(p1 | beta);
  target += dirichlet_lpdf(p2 | beta);
  target += multinomial_lpmf(x | q);
}
