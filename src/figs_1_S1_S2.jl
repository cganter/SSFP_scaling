#= 
==================================================================  
Script, used to generate figures 1, S1, S2

It can be modified to checkout different settings.
(like different choices for R1, R2, n, ϑ, ϕ_inc_deg)
==================================================================  
=#

using SSFP_scaling
using Plots: savefig, pyplot, gr

#-----------------------------------------------------------------
# Choose the backend
#-----------------------------------------------------------------

backend = :gr
# backend = :pyplot

backend == :gr && gr()           # for quick results (also shows categorical heatmaps, which I could not get to work in PyPlot..)
backend == :pyplot && pyplot()   # slow, but somewhat nicer (used for the publication)

#-----------------------------------------------------------------
# Choose the figure
#-----------------------------------------------------------------

make_fig = false # store the figure

fig = :1      # Ernst formula and bSSFP
# fig = :S1     # various coherent configurations
# fig = :S2     # RF spoiling  

#-----------------------------------------------------------------
# Sequences to simulate
#-----------------------------------------------------------------

seqs = []

if fig == :1

    push!(seqs, Ernst())

    # For bSSFP, ϑ = π corresponds to the center of the bSSFP profile.
    # (e.g. spins on resonance with alternating RF pulses)
    # The banding artifact would correspond to ϑ = 0 (or any multiple of 2π)

    push!(seqs, Balanced(Dict(:ϑ => π))) 

elseif fig == :S1

    for n in -2:1
        push!(seqs, Coherent(Dict(:n => n)))
    end

elseif fig == :S2
    
    push!(seqs, Ernst())                               # Ernst formula
    push!(seqs, RFSpoiled(Dict(:ϕ_inc_deg => 117.0)))  # original proposal by Zur et al.
    push!(seqs, RFSpoiled(Dict(:ϕ_inc_deg => 50.0)))   # Siemens
    push!(seqs, RFSpoiled(Dict(:ϕ_inc_deg => 5.0)))    # partial spoiling
    # push!(seqs, RFSpoiled(Dict(:ϕ_inc_deg => 0.0)))    # no spoiling
    # push!(seqs, Coherent(Dict(:n => 0)))               # same as no spoiling

end

#-----------------------------------------------------------------
# Flip angles and relative B1+
#-----------------------------------------------------------------

αs = deg2rad.(range(1, 90, 90))
βs = [0.5, 0.75, 1.0, 1.5, 2]

#-----------------------------------------------------------------
# Generate the plots
#-----------------------------------------------------------------

show_scaling_αβ(seqs, αs, βs, R1 = 0.01, R2 = 0.1)

#-----------------------------------------------------------------
# Save them, if desired (eps format only supported by PyPlot backend)
#-----------------------------------------------------------------

if make_fig && backend == :pyplot
    fig == :1 && savefig("figure_1.eps")
    fig == :S1 && savefig("figure_S1.eps")
    fig == :S2 && savefig("figure_S2.eps")
end
