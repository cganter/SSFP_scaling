#= 
==================================================================  
Testing of partial derivatives from module SSFP_scaling

This will be done graphically, by displaying the deviation
between the exact solution and a linear expansion with respect
to any degree of freedom:

g(h) = |[f(x+h) - f(x)]/h - f'(x)|

f'(x) is correct, if the g(h) vanishes at least linearly.
(∝ h^(n-1), if n denotes the smallest nonzero derivative order ≥ 2)

A different test is needed, if f(x) depends only linearly on x
(like for the proton density prefactor ρ), 
since all higher order derivatives vanish.
In this case, we simply confirm g(h) ≈ 0 without plotting anything.
==================================================================  
=#

using SSFP_scaling
using Random
using Plots: plot, scatter!, gr
gr()

# For the article, we use a fixed seed to remain reproducible.
# To checkout different parameters, just set the following value to false.

fixed_seed = false

rng = fixed_seed ? MersenneTwister(42) : MersenneTwister()

# set the parameter combinations to be simulated

N = 100
α_nom = rand(rng, deg2rad.(range(10, 90, N)))
TR = 1.0
ρ = rand(rng, range(0.5, 1.5, N))
β = rand(rng, range(0.5, 1.5, N))
R1 = rand(rng, range(0.001, 1, N))
R2 = rand(rng, range(R1, 1, N))
ϑ = rand(rng, deg2rad.(range(0, 360, N)))
ϕ_inc_deg = rand(rng, (117.0, 50.0, 150.0))  # Zur et al., Siemens, Philips

# define the invariant part of the dictionary

d = Dict{Symbol,Any}(
        :α_nom => α_nom,
        :TR => TR,
        :ρ => ρ,
        :β => β,
        :R1 => R1,
        :R2 => R2)

# Derivatives are evaluated for the magnetization densities of all SSFP variants.
# (This automatically also tests the underlying derivatives of the various auxiliary functions.)

seqs = []

push!(seqs, Ernst(d))

for n in -3:2
    push!(seqs, par!(Coherent(d), :n, n))
end

push!(seqs, par!(Balanced(d), :ϑ, ϑ))

push!(seqs, par!(RFSpoiled(d), :ϕ_inc_deg, ϕ_inc_deg))

# We test all derivatives, relevant for the CRB.

# Remarks:
# - The Ernst formula does not depend on R2 - i.e. we test, whether zero is returned.
# - The linear dependence on the proton density ρ requires a different approach (cf. introductory remarks).
# - We test Balanced, despite not including it in the CRB. (for reasons, explained in the supplementary information)

der = Dict(
    Ernst => (:ρ, :β, :R1), 
    Coherent => (:ρ, :β, :R1, :R2),
    Balanced => (:ρ, :β, :R1, :R2),
    RFSpoiled => (:ρ, :β, :R1, :R2))

# Now we perform the actual comparisons.
# Note that we make a plot for each s ∈ seqs.

for s in seqs

    typeof(s) == Ernst && @test ∂f(:m,:R2,s) == 0

    plts = [] 
    
    for ∂x in der[typeof(s)]
        
        x = f(∂x,s)
        h = x .* 10 .^range(-5, -1, 10)         # display g(h) logarithmically over 4 decades

        f_appr = f(:m,s) .+ ∂f(:m,∂x,s) .* h    # linear approximation
        par!(s, ∂x, x .+ h)                     # set x + h (which is a vector)
        f_exact = f(:m,s)                       # calculate f(x + h) 
        par!(s, ∂x, x)                          # reset to x (scalar)

        if ∂x == :ρ
            @test f_exact ≈ f_appr              # the linear expansion should hold exactly
        else
            g = abs.((f_exact .- f_appr) ./ h)
            ∂g_∂h = g[1] / h[1]                 # if the first order derivative of g is nonzero, we may estimate the slope
        
            push!(plts, plot(h, ∂g_∂h .* h, xaxis = :log, yaxis = :log, title = string("∂m/∂",∂x), label = "approx"))
            scatter!(h, g, xaxis = :log, yaxis = :log, label = "exact", legend = :topleft)
        end
    end
    
    tit = string(typeof(s))
    typeof(s) == Coherent && (tit = string(tit, ", n = ", f(:n,s)))
    typeof(s) == RFSpoiled && (tit = string(tit, ", ϕ = ", f(:ϕ_inc_deg, s), "°"))

    display(plot(plts..., plot_title = tit))
end