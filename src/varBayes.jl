immutable Options
    maxits::Integer
end

immutable PriorPars
    sigma2_beta::Float64
    a_sigma2::Float64
    b_sigma2::Float64
    a_tau::Array{Float64}
    b_tau::Array{Float64}
    a::Float64
    b::Array{Float64}
end

varBayes = function (dat::Dat, priors::PriorPars, options::Options)
	# load the dat
	y = dat.y
	X = dat.X
	K = dat.K
	Z = dat.Z
	z_P = dat.z_P
	# load the priors
	maxIter = options.maxits
	# load the priors
	sigma2_beta = priors.sigma2_beta
	a_sigma2 = priors.a_sigma2
	b_sigma2 = priors.b_sigma2
	a_tau = priors.a_tau
	b_tau = priors.b_tau
	a = priors.a
	b = priors.b

	# Estimation for Beta
	beta_E = rand(1:10, 2)
	# Estimation for sigma2
	sigma2_E = b_sigma2/(a_sigma2 - 1)
	# Estimation for tau
	tau_E = fill(2.1, K) #b_tau/(a_tau - 1)
	tau_ln_E = fill(0, K)
	# Estimation for Kernels and keep them fixed
	inv_Sigma = []
	for k in 1:K
	    jgl_results = inv(cov(Z[:, find(z_P .== k)]))
	    inv_Sigma  = push!(inv_Sigma, jgl_results)

	end
	Ks = zeros(K, n, n)
	for k in 1:K
	    Ks[k, :, :] = calcSigma(Z[:, find(z_P .== k)], inv_Sigma[k])
	end
	# Estimation for H
	H_E = zeros(n, K);
	H_V = zeros(K, n, n);
	for k in 1:K
	    sigma_k = tau_E[k] * Ks[k, :, :]
	    # add a diagonal matrix to the original matrix to make sigma p.d.
	    sigma_k = max(0, -minimum(eigvals(sigma_k))) * I
	    #sigma_k = cholfact(Hermitian(sigma_k))
	    h = rand(MvNormal(tau_E[k] * Ks[k, :, :]), 1)
	    H_E[:, k] = h
	    H_V[k, :, :] = eye(n)
	end
	# Estimation for gamma
	gamma_E = fill(0.5, K)

	"variational Bayes starts now"
	for iter in 1:maxIter
	    # Report progress of mcmc iter
	    if (mod(iter, 500) == 0)
	        println("iter = $iter\n")
	    end

	    # Update Beta
	    beta_V = inv(1/sigma2_beta * eye(p_beta) + 1/sigma2_E * X' * X)
	    beta_E = 1/sigma2_E * beta_V * X' * (y - sum(H_E, 2))
	    #println("beta_E = $beta_E\n")

	    # Update sigma2
	    sum_h = 0
	    for k in 1:K
	        sum_h = sum_h + H_E[:, k]' * H_E[:, k] +  trace(H_V[k, :, :]) - 2 * H_E[:, k]' * y;
	    end
	    idx =  collect(combinations(1:K, 2))
	    for j in 1:length(idx)
	        sum_h = sum_h + 2 * H_E[:, idx[j][1]]' * H_E[:, idx[j][2]]
	    end
	    sum_h = sum_h +  y' * y
	    A = (X * beta_E)' * (X * beta_E) + trace(X * beta_V * X') - 2 * beta_E' * X' * (y - sum(H_E, 2)) + sum_h
	    shape_sigma2 = n/2 + a_sigma2
	    rate_sigma2 = 1/2 * A + b_sigma2
	    sigma2_E = (rate_sigma2/(shape_sigma2 - 1))[1, 1];
	    #println("sigma2_E = $sigma2_E\n")

	    for k in 1:K
	        # Update gamma
	        num  = - 1/2 * 1/sigma2_E * A
	        col_index = setdiff(1:K, k)
	        sum_h = 0
	        for j in col_index
	            sum_h = sum_h + H_E[:, j]' * H_E[:, j] + trace(H_V[j, :, :]) - 2 * H_E[:, j]' * y
	        end
	        idx =  collect(combinations(1:K, 2))
	        for j in 1:length(idx)
	            sum_h = sum_h + 2 * H_E[:, idx[j][1]]' * H_E[:, idx[j][2]]
	        end
	        sum_h = sum_h + y' * y
	        B =  (X * beta_E)' * (X * beta_E) + trace(X * beta_V * X') - 2 * beta_E' * X' * (y - sum(H_E, 2)) + sum_h
	        den = - 1/2 * 1/sigma2_E * B
	        lik_odd = exp(num - den)
	        sigm_1 = exp(a + sum(b[k, col_index] .* gamma_E[col_index]))/ (1 + sum(b[k, col_index] .* gamma_E[col_index]))
	        sigm_0 = 1/ (1 + exp(a + sum(b[k, col_index] .* gamma_E[col_index])))
	        prior_odd = sigm_1/sigm_0
	        odd = (lik_odd * prior_odd)[1,1]
	        gamma_E[k] = odd/(odd + 1)

	        # Update tau
	        shape_tau = n/2 + a_tau[k]
	        B = H_E[:, k]' * inv(Ks[k, :, :]) * H_E[:, k]
	        + trace(inv(Ks[k, :, :]) * H_V[k, :, :])
	        rate_tau = 1/2 * B[1, 1] + b_tau[k]
	        tau_E[k] = rate_tau/(shape_tau - 1)
	        tau_E[k] = gamma_E[k] * tau_E[k] + 0

	        # Update H
	        col_index = setdiff(1:K, k)
	        V_h = inv((1/sigma2_E) * eye(n) + inv(tau_E[k] * Ks[k, :, :]))
	        m_h = 1/sigma2_E * V_h * (y - X * beta_E - sum(H_E[:, col_index], 2))
	        H_E[:, k] = gamma_E[k] * m_h + 0
	        H_V[k, :, :] = gamma_E[k] * gamma_E[k] * V_h + 0
	    end
	    #println("gamma_E = $gamma_E\n")
	    #println("tau_E = $tau_E\n")
	    #println("\n")
	end

	return beta_E, sigma2_E, gamma_E, H_E
end
