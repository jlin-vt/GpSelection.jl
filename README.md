# GpSelection

[![Build Status](https://travis-ci.org/jlin-vt/GpSelection.jl.svg?branch=master)](https://travis-ci.org/jlin-vt/GpSelection.jl)

[![codecov.io](http://codecov.io/github/jlin-vt/GpSelection.jl/coverage.svg?branch=master)](http://codecov.io/github/jlin-vt/GpSelection.jl?branch=master)

# GpSelection
An implementation of Gaussian process selections in semiparametric regression for multi-pathway analysis in Julia.

## Installation

```
Pkg.clone("git@github.com:jlin-vt/GpSelection.jl.git")
Pkg.build("GpSelection")
```

## An Example
```
n = 100;  K = 5; num_sig = 2; df = 1.0; p_beta = 2
dat = generateData(K, p_beta, n, df, num_sig)
options = Options(100)
priors = PriorPars(5, 10, 10, fill(5, K), fill(10, K), 0.0, zeros(K, K))

@elapsed beta_E, sigma2_E, gamma_E, H_E = varBayes(dat, priors, options)

```

## To Do
* Allow user to define the location of significant pathways instead of putting them at the end.
* More coming up on unit testing.
