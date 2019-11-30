type Z_data
    data::Array{Float64}
    Omega::Array{Float64}
end

type Dat
    y::Array{Float64}
    X::Array{Float64}
    K::Integer
    Z::Array{Float64}
    z_P::Array{Integer}
    sigma2_true::Float64
    beta_true::Array{Float64}
    gamma_true::Array{Integer}
    tau_true::Array{Float64}
end

"""
    generateZ(p, n, signal)

Generate geneteic data by Gaussian graphical model.

# Arguments
- `p::Integer`: dimension of each pathway.
- `n::Integer`: sample size.
- `signal::Bool`: if to generate functional pathway.

# Examples
```julia-repl
julia> Z =  generateZ(p, n, signal)
```
"""
function generateZ(p::Integer, n::Integer, signal::Bool)
	if signal
        Omega = matrixdepot("toeplitz", [1; .5; .25; fill(0, p - 3)]);
        # make Omega invertible
        Omega = fixMatrix(Omega, 2.0)
        @assert all(eig(Omega)[1] .> 0)
        # generate data
        Cov_True = inv(Omega);
        # make it symmetric
        for i = 1:(size(Cov_True, 1) - 1)
            for j = (i + 1): size(Cov_True, 1)
                Cov_True[j, i] = Cov_True[i, j]
            end
        end
        data = rand(MvNormal(Cov_True), n)'
	else
        Cov_True = eye(p) * 1/10000
        data = rand(MvNormal(Cov_True), n)'
        Omega = inv(Cov_True)
	end
	return Z_data(data, Omega)
end

"""
    generateData(p, n, signal)

Generate geneteic data and clinical data.

# Arguments
- `K::Integer`: dimension of each pathway.
- `p_beta::Integer`: number of clinical covariates.
- `n::Bool`: sample size.
- `df::Integer`: scale parameter.
- `num_sig::Bool`: number of significant pathways to include.

# Examples
```julia-repl
julia> generateData(K, p_beta, n, df, num_sig)
```
"""

function generateData(K::Integer, p_beta::Integer, n::Integer, df::Float64, num_sig::Integer)
    # set the true parameters
    sigma2_true = 1.0;
    p_tau = fill(20, K);
    beta_true = rand(Uniform(0, 10), p_beta);
    gamma_true = [fill(0, K - num_sig); fill(1, num_sig)]
    tau_true = gamma_true * df
    signal_true = gamma_true .== 1
    # generate Z(list)
    Z = []
    for k in 1:K
        Z = push!(Z, generateZ(p_tau[k], n, signal_true[k]))
    end
    # X is correlated with first 2 gens in the last pathway
    X = 3 * cos(Z[K].data[:, 1:p_beta]) .+ rand(Normal(), n)
    # generate random effects h
    H_true = zeros(n, K)
    for k in 1:K
        if gamma_true[k] != 0
            H_true[:, k] = rand(MvNormal(tau_true[k] * calcSigma(Z[k].data, Z[k].Omega)), 1)
        end
    end
    # give pathway membership and overwrite Z with z
    z_P = Int64[]
    z = zeros(n, sum(p_tau))
    for k in 1:K
        z_P = append!(z_P, fill(k, p_tau[k]))
        z[:, find(z_P .== k)] =  Z[k].data
    end
    # combine fixed and random effects
    y = X * beta_true + sum(H_true, 2) + rand(Normal(), n)
	return Dat(y, X, K, z, z_P, sigma2_true, beta_true, gamma_true, tau_true)
end
