module fdeSolver

using SpecialFunctions

greet() = print("Hey, let's solve some FDEs!")

include("../FDE_solver/improveit10.jl")
include("../FDE_solver/defineY.jl")
include("../FDE_solver/a_n0_function.jl")
include("../FDE_solver/alpha_function.jl")
include("../FDE_solver/phi_function.jl")
include("../FDE_solver/gamma_function.jl")

export(improveit10)

end
