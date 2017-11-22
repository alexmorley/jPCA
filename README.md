## READ THIS:
** Work in Progress. I'm almost 100% sure this isn't correct yet! You have been warned! **

# jPCA
Julia implementation of jPCA - a PCA variant that captures rotational dynamics

Original methods are shown in [Neural population dynamics during reaching. Churchland ... Shenoy 2012](https://images.nature.com/original/nature-assets/nature/journal/v487/n7405/extref/nature11129-s1.pdf)

and matlab code at [Churchland Lab Site](http://churchlandlab.neuroscience.columbia.edu/links.html)

# Requirements
- Should fit into JuliaStats Framework. i.e. main call should be `fit(JPCA::Type, X)

# Steps:
1. PCA (optional): reduce the dimensionality of the space so as not to capture rotational dynamics with very low varaince.
2. Get the difference between all time points `XD = X[T,:] .- X[T.-1,:]`
3. Find regression between previous and difference between points in `X` `M = XD/X[T.-1,:]`
4. Now we want to find the same regression but constrain `Mskew` to be skew-symmetric (e.g. opposite triangle is negative transpose of other and diagonal is zero  or `Mskew = -Mskew'`).
5. The eigenvectors from the decompositon of Mskew can then be recombined to find the planes of this projection, which we can then use in the same way as PCA to project the data from the higher-D space on this lower-D rotation of the PCA space.
6. "The magnitude of the eigenvalues allow us to select the plane (or planes) with the highest frequency and
consistency in the rotational linear dynamical system"
