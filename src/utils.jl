"Fix to ensure positive definiteness by dividing each off-diagonal element by sum of absolute values of off-diagonal elements in its row"
function fixMatrix(A::Array{Float64}, denom_factor::Float64)
	p = size(A, 1);
	for cur_row in 1:p
		cur_sum = sum(abs(A[cur_row, :])) - 1
		if cur_sum != 1
			A[cur_row, :] = A[cur_row, :] / (denom_factor * cur_sum);
		end
		# Make sure diagonal entries are still 1
		A[cur_row, cur_row] = 1;
	end
	# Final matrix is average of matrix with its transpose
	A = (A + A')/2
    return A
end

"Calculate Covariance Matrix using Gaussian Kernel"
function calcSigma(X, inv_Sigma)
	Sigma = zeros(size(X, 1), size(X, 1))
	for i in 1:size(Sigma, 1)
	  for j in 1:size(Sigma, 2)
          Sigma[i, j] = exp(- 1/2 * (X[i, :] - X[j, :])' * inv_Sigma * (X[i, :] - X[j, :]))[1]
	  end
	end
	return Sigma
end
