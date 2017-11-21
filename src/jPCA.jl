module jPCA

export fit, JPCA

using MultivariateStats,Optim

include("reshape_skew.jl")

import MultivariateStats.fit

type JPCA 
    pca
    M
    Mskew
    proj
end

"""
    fit(jpca::Type{JPCA}, X, k=size(X,2))
NB: X should be a time X features array.
"""
function fit(jpca::Type{JPCA}, X, k=size(X,2))
    k > size(X,2) && error("k must be <= the number of columns in X")
    
    pca   = fit(PCA, X, pratio=1., maxoutdim=k)
    X2    = pca.proj[:,1:k]
    
    T     = (1:size(X2,1))[2:end]
    Tpre  = (1:size(X2,1))[1:end-1]
    dX2 = X2[T,:].-X2[Tpre,:]
    M = dX2/X2[Tpre,:]

    Mskew = skewsymregress(dX2, X2[Tpre,:]) 

    JPCA(pca, M, Mskew, NaN)
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
