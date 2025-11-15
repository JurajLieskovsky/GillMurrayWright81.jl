module GMW81

using LinearAlgebra

"""
Computes the positive-definite factorizaton of `A`.

Returns the permutation `p` and lower factor `L` where `A[p,p] ≈ L * L'` if `A` is positive-definite.
If `A` is not positive definite `(A + E)[p,p] = L * L'` where `E` is a positive diagonal matrix.

The optional argument `δ` is the smallest possible value on the diagonal of `L`. By default, `δ` is
the machine precision of the element type of `A`.

"""
function factorize(A::AbstractMatrix{T}, δ=eps(T)) where {T<:AbstractFloat}
    n, m = size(A)
    @assert n > 1
    @assert n == m

    # result allocation
    p = Vector{Int}(undef, n)
    L = Matrix{T}(undef, n, n)

    # inital copying
    p .= 1:n
    L .= A

    # calculation of β²
    γ = mapreduce(e -> abs(e), max, view(A, diagind(A)))
    ξ = mapreduce(e -> abs(e), max, LowerTriangular(view(A, 2:n, 1:n-1)))
    β2 = max(γ, ξ / sqrt(n^2 - 1), eps(T))

    for j in 1:n
        # find largest diagonal element in the block to be factorized
        q = j
        for i in j:n
            abs(L[i, i]) >= abs(L[q, q]) && (q = i)
        end

        # swap permutation
        p[q], p[j] = p[j], p[q]

        # swap rows
        for i in 1:n
            L[q, i], L[j, i] = L[j, i], L[q, i]
        end

        # swap cols
        for i in j:n
            L[i, q], L[i, j] = L[i, j], L[i, q]
        end

        # calculate factorization
        θ_j = j == n ? 0 : maximum(view(L, j+1:n, j))
        d_j = max(δ, abs(L[j, j]), θ_j^2 / β2)

        L[j, j] = sqrt(d_j)

        if j < n
            view(L, j, j+1:n) .= 0
            view(L, j+1:n, j) .= view(L, j+1:n, j) / L[j, j]
            view(L, j+1:n, j+1:n) .-= view(L, j+1:n, j) * view(L, j+1:n, j)'
        end
    end

    return p, L
end


"""
Reconstruct a positive-definite matrix `A` from the permutation `p` and lower factor `L`.

"""
function reconstruct!(A::AbstractMatrix{S}, p::Vector{Int}, L::Matrix{T}) where {S<:AbstractFloat,T<:AbstractFloat}
    mul!(view(A, p, p), L, L')
    return nothing
end


"""
Reconstruct a positive-definite matrix `A` from the permutation `p` and lower factor `L`.

"""
function reconstruct(p::Vector{Int}, L::Matrix{T}) where {T<:AbstractFloat}
    A = zeros(T, size(L)...)
    reconstruct!(A, p, L)
    return A
end

end # module GMW81

