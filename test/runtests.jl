using Test
using LinearAlgebra
using GillMurrayWrightFactorization: factorize, reconstruct, reconstruct!

# Factorization of a positive-definite matrix
begin
    A = [
        1.0 0.3 0.0
        0.3 2.0 -0.2
        0.0 -0.2 3.0
    ]

    F = factorize(A)

    # test that the reconstructed matrix
    # is the same as the original

    ## out-of-place reconstruction
    @test isapprox(A, reconstruct(F))

    ## in-place reconstruction
    reA = zeros(Float64, 3, 3)
    reconstruct!(reA, F)

    @test isapprox(A, reA)
end

# Factorizaton of an indefinite matrix
begin
    A = [
        1.0 0.3 0.0
        0.3 2.0 -0.2
        0.0 -0.2 -3.0
    ]

    # test that the original matrix is not posdef
    @test !isposdef(A)

    # factorize and reconstruct matrix
    F = factorize(A)
    reA = reconstruct(F)

    # test that the reconstructed matrix is posdef
    @test isposdef(reA)

    # test that only diagonal entries were modified
    A[diagind(A)] .= 0
    reA[diagind(reA)] .= 0

    @test isapprox(A, reA)
end
