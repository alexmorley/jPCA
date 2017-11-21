module jPCA

export fit, JPCA

using MultivariateStats,Optim

include("reshape_skew.jl")
include("getRealVs.jl")

import MultivariateStats.fit

type JPCA 
    pca
    M
    Mskew
    proj
end

"""
    fit(jpca::Type{JPCA}, X, k=10 or num features)
NB: X should be a time X features array.
"""
function fit(jpca::Type{JPCA}, X, k=size(X,2) < 10 ? size(X,2) : 10)
    size(X,2) > size(X,1) && error("Must have more rows than columns")
    k > size(X,2) && error("k must be <= the number of columns in X")
	isodd(k) && error("k must be even")   
 
    pca   = fit(PCA, X, pratio=1., maxoutdim=k)
    X2    = pca.proj[:,1:k]
    
    T     = (1:size(X2,1))[2:end]
    Tpre  = (1:size(X2,1))[1:end-1]
    dX2 = X2[T,:].-X2[Tpre,:]
    M = dX2/X2[Tpre,:]

    Mskew = skewsymregress(dX2, X2[Tpre,:]) 
    D,V   = eig(Mskew)
    V     .= V[:,sortperm(abs.(D))]

    jPCs = zeros(size(V)...)
	for pair in 1:3
		vi1 = 1+2(pair-1)
		vi2 = 2pair

		VconjPair  = V[:, [vi1,vi2]]
		evConjPair = D[[vi1, vi2]]
		jPCs[:,vi1] .= getRealVs(VconjPair, evConjPair)[1]
		jPCs[:,vi2] .= getRealVs(VconjPair, evConjPair)[2]
    end
	proj = pca.proj * jPCs 

    JPCA(pca, M, Mskew, proj)
end

function skewsymregress(dX,X)
    M0 = X\dX
    M0 .= 0.5*(M0 - M0')
    m0 = reshape_skew(M0)

    function skewsymderiv(m)
        # Evaluate objective function
        f = vecnorm( dX - X*reshape_skew(m))^2;

        # Evaluate derivative
        # D = (dX - X*reshape_skew(m))'*X;
        # df = 2*reshape_skew(D - D')
        f
    end

    M = optimize(skewsymderiv, m0)
    return reshape_skew(M.minimizer)
end

end
