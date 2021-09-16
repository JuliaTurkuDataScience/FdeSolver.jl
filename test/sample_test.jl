using Test
using FdeSolver

@testset "FdeSolver.jl" begin
    b = FdeSolver.greet()
    @test @isdefined(b)
end
