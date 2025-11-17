# GillMurrayWright81.jl

This package provides a basic implementation of the `GMW81` algorithm, described in [Practical Optimization](https://doi.org/10.1137/1.9781611975604) by P. E. Gill, W. Murrray and M. H. Wright. The implementation itself closely follows the later description of the algorithm in [What to Do When Your Hessian is Not Invertible](https://doi.org/10.1177/0049124103262681) by J. Gill and G. King.

Given a symmetric matrix `A` the algorithm performs a modified Cholesky decomposition (including pivoting) of `A+E`, where `E` is a positive diagonal matrix, designed as the minimal necessary correction required to make `A+E` positive-definite.

Anecdotally, this standard algorithm has significantly better convergence properties on certain problems when compared to the approach taken in [PostiveFactorizations.jl](https://github.com/timholy/PositiveFactorizations.jl). The option to compare the two was the primary motivation behind the creation of this package. As an added benefit, one may simply cite the first reference when using the package for academic research.

## Usage

A symmetric matrix `A::AbstractMatrix{T}` where `T<:AbstractFloat` can be factorized using
```
F = GillMurrayWright81.factorize(A)
```
or
```
F = GillMurrayWright81.factorize(A, δ) 
```
where `p::Vector{Int}` is the permutation and `L::Matrix{T}` the lower triangular factor. The optional argument `δ` (which defaults to `eps(T)`) is the smallest value permitted on the diagonal of `L`. To reconstruct a positive-definite matrix `Ã` you may use either
```
Ã = GillMurrayWright81.reconstruct(p, L)
```
or the in-place version
```
GillMurrayWright81.reconstruct!(Ã, p, L)
```
It is important to note that `factorize` only uses the lower triangular elements of `A` to compute the positive-definite factorization and does not check if `A` is symmetric.
