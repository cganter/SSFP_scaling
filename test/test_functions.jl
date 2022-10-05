#= 
==================================================================  
Testing of functions from module SSFP_scaling
==================================================================  
=#

using SSFP_scaling

# check types

@test Ernst <: SSFP
@test Configuration <: SSFP
@test Coherent <: Configuration
@test RFSpoiled <: Configuration
@test Balanced <: SSFP

# define a few parameters for the simulations

ρ = 0.8
β = 0.9
α_nom = collect(range(1, π/2, 90))
α = β .* α_nom
R1, TR = 0.001, 1.0

# test par! (set and overwrite parameters)

is = Ernst()
par!(is, :TR, 2TR)
par!(is, Dict(:R1 => 0.1R1, :α_nom => 2α_nom, :β => β))
par!(is, :α_nom, α_nom)
par!(is, Dict(:R1 => R1, :TR => TR, :ρ => ρ))

@test f(:TR,is) == TR
@test f(:R1,is) == R1
@test f(:α_nom,is) == α_nom
@test f(:β,is) == β
@test f(:ρ,is) == ρ

# calculate a and b, check their correct values and that they are 
# stored in the SSFP instance

E1 = exp(- R1 .* TR)
cα = cos.(α)
sα = sin.(α)

a = f(:a,is)
b = f(:b,is)

@test a ≈ 1.0 .- E1
@test :a ∈ keys(par(is))
@test b ≈ 1.0 .- cα
@test :b ∈ keys(par(is))

# resetting basic parameters invalidates derived ones (in this case a and b)

par!(is, :R1, R1)

@test :a ∉ keys(par(is))
@test :b ∉ keys(par(is))

# check, whether the Ernst formula is reproduced correctly.

ErnstFormula = ρ .* sα .* (1.0 .- E1) ./ (1.0 .- cα .* E1)

@test f(:m,is) ≈ ErnstFormula

# Coherent configurations are tested indirectly:
# We make use of the fact that balanced SSFP (bSSFP) can be expressed in terms of all configurations
# and the effective off-resonance ϑ, as explained in the supplementary information.
# Since this is just how the type Balanced is implemented, we compare Balanced against an available
# closed form analytical bSSFP steady state expression.

# Effective precession angle (off-resonance - RF phase increment).
# A whole range is tested in order to check the individual configurations.

ϑ = collect(range(-π, π, 73))'  
cϑ = cos.(ϑ)
sϑ = sin.(ϑ)
eiϑ = exp.(im .* ϑ)

# transverse relaxation

R2 = 0.01
E2 = exp(- R2 .* TR)

# Closed form solution of the bSSFP steady state

D = (1.0 .- E1 .* cα) .* (1.0 .- E2 .* cϑ) .- (E1 .- cα) .* (E2 .- cϑ) .* E2
bSSFP_theo = ρ .* sα .* (1.0 .- E1) .* (1.0 .- E2 .* (cϑ .- im .* sϑ)) ./ D

# FID and ECHO configurations

fid = par!(Coherent(par(is)), Dict(:R2 => R2, :n => 0))
echo = par!(Coherent(par(fid)), :n, -1)

# Explicit expression in terms of FID, ECHO, A and off-resonance (ϑ).

bSSFP_expl = f(:m,fid) ./ (1.0 .- f(:A,fid) .* eiϑ) .+ f(:m,echo) ./ (eiϑ .- f(:A,fid))

# bSSFP steady state as returned (and to be tested) in SSFP_scaling.

bSSFP = f(:m, par!(Balanced(par(fid)), :ϑ, ϑ))

# All three expression must match (approximately).

@test bSSFP_theo ≈ bSSFP_expl ≈ bSSFP

# Strictly speaking, we just confirmed the cases FID (n == 0) and ECHO (n == -1) by now.
# Also, we know that the expression for A seems to be correct.
# Therefore, we now confirm higher order configurations explicitly in terms of FID, ECHO and A.

cfg = Coherent(par(fid))

for n in 1:10
    @test f(:m, par!(cfg, :n, n)) ≈ f(:m,fid) .* f(:A,fid).^n
    @test f(:m, par!(cfg, :n, -(n+1))) ≈ f(:m,echo) .* f(:A,fid).^n
end
