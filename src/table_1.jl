#= 
==================================================================  
Script, used to generate table 1

For given sets of steady state SSFP acquisitions, 
we assess the attainable precision for estimators (R1, R2), 
depending on whether the relative B1+ (β) is known or not.

(In the latter case, we assume β to be fitted from the available data as well.)

The script allows for a flexible setup of acquisitions:

    - number of spoiled and coherent scans with respective
    - free choice between ideal spoiling (Ernst formula) and actual RF spoiling
    - arbitrary flip angles
    - arbitrary configuration orders (coherent only)

and other relevant parameters (R1, R2, ...).

Finally, the script allows to validate the calculated CRB results
via a Monte Carlo (MC) simulation:

To do so, it suffices to provide a nonzero value for the number (n_MC)
of MC samples below. The agreement between CRB and MC improves with 
smaller σ (like 1e-4), i.e. higher SNR, since the quadratic expansion becomes 
less reliable with increasing deviation from the minimum. 

The chosen noise value σ = 1e-3 roughly corresponds to a SNR of about 100, since 
the SSFP signal is of order 0.1 in units of the proton density ρ.

**** Note: *****
The MC simulation is rather time consuming and therefore disabled by default.
The time-critical part exploits multi-threading though.
To benefit from this, just start julia with a suitable number of threads:

julia -t <number of threads>

or simply (if supported)

julia -t auto
==================================================================  
=#

using SSFP_scaling
using LinearAlgebra, Printf, Random, Statistics, Optim, NLSolversBase, ProgressMeter

#-----------------------------------------------------------------
# number of MC samples
#-----------------------------------------------------------------

n_MC = 0  # no MC simulations
# n_MC = 1000 (used for table 1)

#-----------------------------------------------------------------
# Fixed parameters (change as desired)
#-----------------------------------------------------------------

d = Dict(:ρ => 1.0, :R1 => 0.01, :R2 => 0.1, :β => 1.0, :TR => 1.0)

#-----------------------------------------------------------------
# Ernst angle for the given R1 and TR
#-----------------------------------------------------------------

α_Ernst = acos(exp(- d[:R1] * d[:TR]))    # Ernst angle

#-----------------------------------------------------------------
# data noise [ρ] * ρ (assumed to be the same for each acquisition)
#-----------------------------------------------------------------

σ = 1e-3 .* d[:ρ]

#-----------------------------------------------------------------
# RF spoiled acquisitions
# (Check out different settings by adjusting α_rel)
#-----------------------------------------------------------------

# α_rel = []            # no (RF) spoiled acquisition
α_rel = [1.0]         # one measurement (at the Ernst angle)
# α_rel = [0.5, 1.0]    # two measurements (half and full Ernst angle)

α_sp = α_rel .* α_Ernst           # flip angles for spoiled acquisitions

#-----------------------------------------------------------------
# Coherent configurations
# (Check out different settings of α_co and n_co)
#-----------------------------------------------------------------

# multiple flip angles, FID only

α_co = deg2rad.([10.0, 30.0, 50.0])
n_co = [0, 0, 0]

# constant flip angle, various configuration orders (similar to TESS)
# 
# ***** Note: *****
# If β is considered unknown, this variant requires at least two (RF) spoiled acqisitions.
# Otherwise the rank of the matrix ∂m_∂θ (and thus the Fisher matrix)
# would be smaller than length(θ).
# This can be proven by making use of the following facts (cf. SI):
# - m^(n) is a function of c, d, mE.
# - The latter are constant for constant flip angles

# α_co = deg2rad.([30.0, 30.0, 30.0])
# n_co = [-1, 0, 1]

#-----------------------------------------------------------------
# Apply settings and set up the scans
#-----------------------------------------------------------------

seqs = []    # sequence vector

for α in α_sp   # (RF) spoiled SSFP

    # ideal spoiling (Ernst formula)
    push!(seqs, par!(Ernst(d), Dict(:α_nom => α)))

    # real RF spoiling (should yield similar results)
    # push!(seqs, par!(RFSpoiled(d), Dict(:α_nom => α, :ϕ_inc_deg => 50.0)))

end

for (n, α) in zip(n_co, α_co)   # coherent, unbalanced SSFP

    push!(seqs, par!(Coherent(d), Dict(:α_nom => α, :n => n)))

end

#-----------------------------------------------------------------
# We estimate the Cramer Rao bound (CRB) for two cases:

# - β is assumed to be known and set correctly from a separate scan
# - β is not known and fitted from the SSFP scans
#-----------------------------------------------------------------

θs_β_known = (:R1, :R2, :ρ)
∂m_∂θ_β_known = ∂m_∂θ(seqs, θs_β_known)
rank_ok_β_known = rank(∂m_∂θ_β_known) >= length(θs_β_known)
rank_ok_β_known && (crb_β_known = crb(seqs, θs_β_known, σ))

θs_β_unknown = (:R1, :R2, :ρ, :β)
∂m_∂θ_β_unknown = ∂m_∂θ(seqs, θs_β_unknown)
rank_ok_β_unknown = rank(∂m_∂θ_β_unknown) >= length(θs_β_unknown)
rank_ok_β_unknown && (crb_β_unknown = crb(seqs, θs_β_unknown, σ))

#-----------------------------------------------------------------
# Output, corresponding to table 1 in the manuscript
#-----------------------------------------------------------------

rank_ok_β_known && (
    relerr_β_known = Dict(
        θ => v for (θ, v) in zip(
            θs_β_known,
            sqrt.(diag(crb_β_known)) ./ [f(θ, seqs[1]) for θ in θs_β_known],
        )
    )
)

rank_ok_β_unknown && (
    relerr_β_unknown = Dict(
        θ => v for (θ, v) in zip(
            θs_β_unknown,
            sqrt.(diag(crb_β_unknown)) ./ [f(θ, seqs[1]) for θ in θs_β_unknown],
        )
    )
)

rank_ok_β_known &&
    rank_ok_β_unknown &&
    (penalty = Dict(θ => relerr_β_unknown[θ] ./ relerr_β_known[θ] for θ in (:R1, :R2)))

@printf "==================\n"
@printf "%d SSFP measurements:\n\n" length(α_sp) + length(α_co)
for (i, α) in enumerate(α_sp)
    @printf "Spoiled FID with α = %.1f°\n" rad2deg(α)
end
for (n, α) in zip(n_co, α_co)
    @printf "Coherent configuration:\t n = %2d, α = %.1f°\n" n rad2deg(α)
end
@printf "\n------------------\n"
@printf "CRB estimates, assuming data noise σ = %.1e * ρ\n\n" σ
@printf "β known and set properly:\n\n"

if rank_ok_β_known
    @printf "rel. error R1 = %.1f %%\n" 100 * relerr_β_known[:R1]
    @printf "rel. error R2 = %.1f %%\n" 100 * relerr_β_known[:R2]
else
    @printf "Fit impossible due to insufficient rank of Fisher matrix\n"
end

@printf "\nβ unknown and fitted:\n\n"

if rank_ok_β_unknown
    @printf "rel. error R1 = %.1f %%\n" 100 * relerr_β_unknown[:R1]
    @printf "rel. error R2 = %.1f %%\n" 100 * relerr_β_unknown[:R2]
else
    @printf "Fit impossible due to insufficient rank of Fisher matrix\n"
end

@printf "\nPrecision penalty:\terror(β fitted) / error(β set properly):\n\n"

if rank_ok_β_known && rank_ok_β_unknown
    @printf "R1 penalty = %.2f\n" penalty[:R1]
    @printf "R2 penalty = %.2f\n" penalty[:R2]
else
    @printf "Not computable due to insufficient rank of (at least one) Fisher matrix"
end

#-----------------------------------------------------------------
# Optional MC simulation, to test whether the CRB agrees with 
# the estimators' covariance matrix, obtained from numerical minimization.
# Only executed, if n_MC > 0 and full rank for both optimizations.
#-----------------------------------------------------------------

if n_MC > 0 && rank_ok_β_known && rank_ok_β_unknown

    ssfp(seqs) = [f(:m, s) for s in seqs]
    ∇ssfp(seqs, θks) = [∂f(:m, ∂x, s) for s in seqs, ∂x in θks]

    """
    χ22_fgh!(f, g, h, θs, θks, seqs, data, σ)

    Calculates value, gradients and Hessian for the set of SSFP acquisitions.
    Avoids repetitive calculations.
    """
    function χ22_fgh!(f, g, h, θs, θks, seqs, data, σ)
        # set actual θs in all sequences
        (s -> par!(s, Dict(θk => θ for (θ, θk) in zip(θs, θks)))).(seqs)
        # calculate and return gradient
        g !== nothing && copyto!(g, ∇ssfp(seqs, θks)' * (ssfp(seqs) .- data) ./ σ^2)
        # calculate and return Hessian
        h !== nothing && copyto!(h, ∇ssfp(seqs, θks)' * ∇ssfp(seqs, θks) ./ σ^2)
        # calculate and return function value
        f !== nothing && return ((ssfp(seqs) .- data)' * (ssfp(seqs) .- data)) ./ 2σ^2
    end

    """
    cov_θ_MC(θks, seqs, σ, n_MC; max_iter = 100)

    Monte Carlo (MC) simulation:

    For the set of sequences and tissue parameters (supplied in `seqs`), 
    `n_MC` noisy data samples are generated and subsequently fitted with
    an Interior point Newton solver.

    Input:
        - `θks`: Vector of estimators ∈ {:R1, :R2, :ρ, :β} to be fitted
        - `seqs`: Sequences, including all parameters (including tissue properties)
        - `n_MC`: Number of MC runs (resuls only enter in case of convegence)
        - `max_iter`: max. allowed number of iterations
    Output:
    The function returns a dictionary with the following keys:
        - `:optim`: Vector of the results returned from `optimize()`
        - `:found_θs`: Vector of estimator vectors
        - `:did_not_converge`: Number of MC runs, which did not converge
        - `:cov`: Covariance matrix of estimators
    """
    function cov_θ_MC(θks, seqs, σ, n_MC; max_iter = 100)
        res = Dict()
        res[:optim] = []
        res[:found_θs] = Array{Float64}(undef, 0, length(θks))
        res[:did_not_converge] = 0

        unperturbed_data = ssfp(seqs)
        initial_θs = [f(θk, first(seqs)) for θk in θks]
        lθ = fill(0, length(θks))
        uθ = fill(Inf, length(θks))

        # We use a fixed seed, to work with equivalent noisy data for
        # β known and unknown. 
        # As a side effect, we reproducible with respect to the article.
        # If this is not desired, just set the following value to false.

        fixed_seed = true

        rng = fixed_seed ? MersenneTwister(42) : MersenneTwister()

        # To allow for multiple threads, we make a few copies

        seqs_vec = [deepcopy(seqs) for _ in 1:n_MC]
        data_vec = [unperturbed_data + σ .* randn(rng, length(ssfp(seqs))) for _ in 1:n_MC]

        df_vec = [TwiceDifferentiable(
            Optim.only_fgh!((F, G, H, x) -> χ22_fgh!(F, G, H, x, θks, seqs_, data_, σ)), initial_θs) 
            for (seqs_,data_) in zip(seqs_vec,data_vec)]
            
        dfc_vec = [TwiceDifferentiableConstraints(lθ, uθ) for _ in 1:n_MC]

        # set up the multi-threaded progress bar

        p = Progress(n_MC, 5)
        update!(p, 0)
           
        # to avoid data races

        l = Threads.SpinLock()
        ii = Threads.Atomic{Int}(0)

        # execute the loop

        Threads.@threads for i in 1:n_MC
            res_optim = Optim.optimize(
                df_vec[i],
                dfc_vec[i],
                initial_θs,
                IPNewton(),
                Optim.Options(iterations = max_iter),
            )

            Threads.atomic_add!(ii,1)

            Threads.lock(l)
            
            if Optim.converged(res_optim)
                push!(res[:optim], res_optim)
                res[:found_θs] = [res[:found_θs]; Optim.minimizer(res[:optim][end])']
            else
                res[:did_not_converge] += 1
            end
            
            update!(p, ii[])
            
            Threads.unlock(l)
        end

        # calculate covariance matrix

        res[:cov] = cov(res[:found_θs])

        return res
    end

    @printf "\n------------------\n"
    @printf "Conducting MC simulation for %d samples\n" n_MC
    @printf "and same data noise as for the CRB: σ = %.1e * ρ\n\n" σ
    @printf "1/2: Assuming that β is known and set properly\n"
    MC_β_known = cov_θ_MC(θs_β_known, seqs, σ, n_MC)
	@printf "%d of %d samples did not converge and were excluded\n\n" MC_β_known[:did_not_converge] n_MC
    @printf "2/2: Assuming that β is unknown and fitted\n"
    MC_β_unknown = cov_θ_MC(θs_β_unknown, seqs, σ, n_MC)
	@printf "%d of %d samples did not converge and were excluded\n\n" MC_β_unknown[:did_not_converge] n_MC

    #-----------------------------------------------------------------
    # Output, corresponding to Table S1 in the manuscript
    #-----------------------------------------------------------------

    relerr_β_known = Dict(
        θ => v for (θ, v) in zip(
            θs_β_known,
            sqrt.(diag(MC_β_known[:cov])) ./ [f(θ, seqs[1]) for θ in θs_β_known],
        )
    )

    relerr_β_unknown = Dict(
        θ => v for (θ, v) in zip(
            θs_β_unknown,
            sqrt.(diag(MC_β_unknown[:cov])) ./ [f(θ, seqs[1]) for θ in θs_β_unknown],
        )
    )

    penalty = Dict(θ => relerr_β_unknown[θ] ./ relerr_β_known[θ] for θ in (:R1, :R2))

    @printf "------------------\n"
    @printf "Results of MC simulation\n\n" 
    @printf "β known and set properly:\n\n"

    @printf "rel. error R1 = %.1f %%\n" 100 * relerr_β_known[:R1]
    @printf "rel. error R2 = %.1f %%\n\n" 100 * relerr_β_known[:R2]

    @printf "β unknown and fitted:\n\n"

    @printf "rel. error R1 = %.1f %%\n" 100 * relerr_β_unknown[:R1]
    @printf "rel. error R2 = %.1f %%\n\n" 100 * relerr_β_unknown[:R2]

    @printf "Precision penalty:\terror(β fitted) / error(β set properly):\n\n"

    @printf "R1 penalty = %.2f\n" penalty[:R1]
    @printf "R2 penalty = %.2f\n" penalty[:R2]
end