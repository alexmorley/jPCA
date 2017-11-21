"""
    reshape_skew(x)
Reshapes a n(n-1)/2 vector to a n by n skew symmetric matrix, or vice versa
"""
function reshape_skew(x::Array{Float64,1})
    # first get the size of the appropriate matrix
    # this should be n(n-1)/2 entries.
    # this is the positive root
    nf = (1 + sqrt(1 + 8*length(x)))/2;
    # error check
    if nf != round(nf) #if not an integer
        error("size of the x vector prevents it from being shaped into a skew symmetric matrix.")
    end
    
    n = floor(Int, nf)
    # now make the matrix
    # initialize the return matrix
    Z = zeros(n,n);

    # and the marker index
    indMark = 1;
    for i = 1 : n-1
        # add the elements as appropriate.
        Z[i+1:end,i] = x[indMark:indMark+(n-i)-1];
        # now update the index Marker
        indMark = indMark + (n-i);
    end

    # now add the skew symmetric part
    Z = Z - Z'
    return Z
end
function reshape_skew(x::Array{Float64,2}) 
    # then we are making a vector from a matrix (note that the 
    # standard convention of lower case being a vector and upper case
    # being a matrix is now reversed).

    # first check that everything is appropriately sized and skew
    # symmetric
    size(x,1) == size(x,2) || error("the matrix x is not square")

    # now check for skew symmetry
    if abs(norm(x + x')) > 1e-8
        # this is not skew symmetric.
        error("the matrix x is not skew-symmetric")
    end
    # everything is ok, so take the size
    n = size(x,1);

    # now make the vector Z
    Z = zeros(div(n*(n-1),2))
    indMark = 1;
    for i in 1:(n-1)
        # add the elements into a column vector as appropriate.
        Z[indMark:indMark+(n-i)-1,1] = x[i+1:end,i];
        # now update the index Marker
        indMark = indMark + (n-i);
    end
    return Z
end
