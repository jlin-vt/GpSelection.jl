using Base.Test
using GpSelection

n = 100;  K = 5; num_sig = 2; df = 1.0; p_beta = 2
dat = generateData(K, p_beta, n, df, num_sig)
options = Options(100)
priors = PriorPars(5, 10, 10, fill(5, K), fill(10, K), 0.0, zeros(K, K))

#@elapsed beta_E, sigma2_E, gamma_E, H_E = varBayes(dat, priors, options)
