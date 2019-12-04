# GpSelection
[![Build Status](https://travis-ci.org/jlin-vt/GpSelection.jl.svg?branch=master)](https://travis-ci.org/jlin-vt/GpSelection.jl)
[![codecov.io](http://codecov.io/github/jlin-vt/GpSelection.jl/coverage.svg?branch=master)](http://codecov.io/github/jlin-vt/GpSelection.jl?branch=master)

## Overview
GpSelection.jl is a Julia package for Gaussian process selections in semiparametric regression for multi-pathway analysis.

## Installation
To download the package, use the code below.
```
Pkg.clone("git@github.com:jlin-vt/GpSelection.jl.git")
Pkg.build("GpSelection")
```

## Usage
```Julia
n = 100;  K = 5; num_sig = 2; df = 1.0; p_beta = 2
dat = generateData(K, p_beta, n, df, num_sig)
options = Options(100)
priors = PriorPars(5, 10, 10, fill(5, K), fill(10, K), 0.0, zeros(K, K))

@elapsed beta_E, sigma2_E, gamma_E, H_E = varBayes(dat, priors, options)
```

## Status
The preprint describing the corncob methodology is available [here](https://github.com/jlin-vt/GpSelection.jl). The manuscript has been submitted to _Journal of Statistical Computation and Simulation_.

## Bug Reports / Change Requests
If you encounter a bug or would like make a change request, please file it as an issue [here](https://github.com/jlin-vt/GpSelection.jl/issues).

## License
The GpSelection.jl package is licensed under the MIT "Expat" License.
