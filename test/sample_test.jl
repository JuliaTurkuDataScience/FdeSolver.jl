using Test
using FdeSolver

@testset "FdeSolver.jl" begin
    a = 2
    b = 3
    @test my_f(a, b) == 5
end

@testset "FdeSolver.jl" begin
    b = FdeSolver.greet()
    @test @isdefined(b)
end
