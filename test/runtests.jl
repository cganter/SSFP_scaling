#= 
==================================================================  
Testing of functions and derivatives from module SSFP_scaling
==================================================================  
=#

using SSFP_scaling
using Test

@testset "SSFP_scaling.jl" begin
    include("test_functions.jl") 
    include("test_derivatives.jl") 
end;
