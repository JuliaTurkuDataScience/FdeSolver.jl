using SafeTestsets

@safetestset "solving 1 FDE" begin include("1FDE_test.jl") end

@safetestset "solving 2 parametric FDEs" begin include("2FDE_par_test.jl") end

@safetestset "solving 3 parametric FDEs with inner function" begin include("3FDE_par_test.jl") end

@safetestset "solving FDE of second order)" begin include("order2_test.jl") end
