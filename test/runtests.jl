using jPCA
using Base.Test

# write your own tests here
X = rand(100,100)
Y = fit(JPCA, X, 10)

@test size(Y.M,1) == size(Y.M,2)

@test_throws ErrorException fit(JPCA, rand(5,5), 6)
