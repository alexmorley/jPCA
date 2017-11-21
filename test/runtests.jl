using jPCA
using Base.Test

# main
X = rand(100,100)
Y = fit(JPCA, X, 10)

@test size(Y.M,1) == size(Y.M,2)

@test_throws ErrorException fit(JPCA, rand(5,5), 6)


# skew sym
import jPCA.reshape_skew
# create a random skew sym matrix
Z = Z = rand(100,100)
Z .= 0.5(Z - Z')
@test Z == reshape_skew(reshape_skew(Z))

@test_throws ErrorException reshape_skew(Z[:,1:99])
Z[1,1] += 0.1
@test_throws ErrorException reshape_skew(Z)
