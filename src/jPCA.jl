module jPCA

export fit, JPCA

using MultivariateStats

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
    
    pca   = fit(PCA, X)
    X2    = pca.proj[:,1:k]
    
    T     = (1:size(X2,1))[2:end]
    Tpre  = (1:size(X2,1))[1:end-1]
    M = (X2[T,:].-X2[Tpre,:])/X2[Tpre,:]

    JPCA(pca, M, NaN, NaN)
end

end
