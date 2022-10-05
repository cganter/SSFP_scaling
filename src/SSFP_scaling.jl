#= 
==================================================================  
Definitions for SSFP functions, used to demonstrate
B1+ scaling and to calculate the Cramer Rao bound
==================================================================  
=#

module SSFP_scaling

using LaTeXStrings
using Plots: plot, plot!, cgrad, heatmap, xlabel!, ylabel!

using Measures

export SSFP, Configuration
export par, par!, remove_derived!, f, ∂f
export ∂m_∂θ, fisher, crb, show_scaling_αβ
export Ernst, Coherent, RFSpoiled, Balanced

#= 
==================================================================  
Abstract SSFP sequences and associated generic functions
==================================================================  
=#

"""
```
SSFP
```
Supertype of any SSFP sequence.
"""
abstract type SSFP end

#-----------------------------------------------------------------

"""
```
Configuration
```
`{Coherent, RFSpoiled} <: Configuration <: SSFP`
"""
abstract type Configuration <: SSFP end

#-----------------------------------------------------------------

"""
```
par(s::SSFP)
```
Return the dictionary of sequence, tissue and derived parameters.

**Note**: Any *concrete* subtype of `SSFP` must contain a field `par`.
"""
par(s::SSFP) = s.par

#-----------------------------------------------------------------

"""
```
par!(s::SSFP, p::Symbol, v)
```
Set/overwrite parameter `p` with value `v`.

**Note**: Applying this function *always* invalidates (i.e. removes) any derived parameters!
# Example
```
par!(s, :R1, 0.01)
```
"""
function par!(s::SSFP, p::Symbol, v)
    par(s)[p] = v
    remove_derived!(s)
end

#-----------------------------------------------------------------

"""
```
par!(s::SSFP, d::Dict)
```
Merge `par(s)` with (a copy of) `d`. (`d` overwrites `par(s)` in case of duplicate keys.) 

**Note**: Applying this function *always* invalidates (i.e. removes) any derived parameters!
# Example
```
par!(s, Dict(:R1 => 0.01, :R2 => 0.1))
```
"""
function par!(s::SSFP, d::Dict)
    merge!(par(s), copy(d))
    remove_derived!(s)
end

#-----------------------------------------------------------------

"""
```
remove_derived!(s::SSFP)
```
Remove any key `∉ (:α_nom, :TR, :ρ, :β, :R1, :R2, :n, :ϑ, :ϕ_inc_deg)` from `par(s)`.
"""
function remove_derived!(s::SSFP)
    (
        k -> (k ∉ (:α_nom, :TR, :ρ, :β, :R1, :R2, :n, :ϑ, :ϕ_inc_deg, :rel_err, :rnd_inc, :max_iter) && delete!(par(s), k))
    ).(keys(par(s)))
    return s
end

#-----------------------------------------------------------------

"""
```
f(k::Symbol, s::SSFP)
```
Return value associated with `k`.

If the key `k` can be found in the dictionary `par(s)`, its value is simply returned.

Otherwise, a function `\$k(s)` is evaluated. The returned result is also stored in `par(s)`
under the key `k`, to avoid repetitive evaluations.
# Example
```
s = Ernst(d)             # d is some dictionary with all required parameters.
a1 = f(:a,s)             # The function a(s) is evaluated, the result stored in par(s)[:a] and returned to a1.
a2 = f(:a,s)             # Since the key :a now exists in par(s), a2 is simply set to par(s)[:a] without calling a(s) again.
par!(s,:R1,0.01)         # The basic parameter :R1 is reset. Since par(s)[:a] is no longer valid, it is removed from s.
a3 = f(:a,s)             # cf. a1
```
"""
f(k::Symbol, s::SSFP) = k ∈ keys(par(s)) ? par(s)[k] : (par(s)[k] = eval(:($k($s))))

#-----------------------------------------------------------------

"""
```
∂f(f_::Symbol, ∂x::Symbol, s::SSFP)
```
Calculates the partial derivative ∂p/∂q(s) and stores the result in `par(s)[:∂p_∂q]`,
given that `f_ == :p` and `∂x == :q`.

For `f_ == :m` the restriction `∂x ∈ (:ρ, :β, :R1, :R2, :a, :b, :E2)` applies.
"""
function ∂f(f_::Symbol, ∂x::Symbol, s::SSFP)
    if f_ == :m && ∂x ∈ (:ρ, :β, :R1, :R2)      # generic cases
        (key_ = Symbol(:∂m_∂, ∂x)) in keys(par(s)) && return par(s)[key_]
        ∂x == :ρ && return (par(s)[key_] = f(:m, s) ./ f(:ρ, s))
        ∂x == :β && return (par(s)[key_] = ∂f(:m, :b, s) .* ∂f(:b, :β, s))
        ∂x == :R1 && return (par(s)[key_] = ∂f(:m, :a, s) .* ∂f(:a, :R1, s))
        ∂x == :R2 && return (par(s)[key_] = ∂f(:m, :E2, s) .* ∂f(:E2, :R2, s))
    else
        ∂f_ = Symbol(:∂, f_)
        (key_ = Symbol(:∂, f_, :_∂, ∂x)) in keys(par(s)) ? par(s)[key_] :
        (par(s)[key_] = eval(:($∂f_($(Meta.quot(∂x)), $s))))
    end
end

#= 
==================================================================  
Testing the validity of the scaling law (1)
==================================================================  
=#

"""
    show_scaling_αβ(seqs, αs, βs;
    ρ = 1.0,
    TR = 1.0,
    R1,
    R2,
    subplt_size = (330,220),
    dpi = 100,
    margin = 2mm,
    percent = true,
    errbnds = (-0.11, 0.11),
    ncats = 11,
    cmp = :vik,
    N = 200)

Testing the validity of the scaling law (1)

We calculate the RHS of Eq. (1) for a range of β.
The calculation for β=1 is also equivalent to the LHS of Eq. (1).
"""
function show_scaling_αβ(seqs, αs, βs;
    ρ = 1.0,
    TR = 1.0,
    R1,
    R2,
    subplt_size = (330,220),
    dpi = 100,
    margin = 2mm,
    percent = true,
    errbnds = (-0.11, 0.11),
    ncats = 11,
    cmp = :vik,
    N = 200)

    αs = reshape(αs, 1, length(αs))
    βs = βs[:]
    βr = range(min(βs...), max(βs...), N)
    errfak = percent ? 100.0 : 1.0

    plts = Matrix(undef, 2, length(seqs))

    for (is,s) in enumerate(seqs)
        if typeof(s) == Ernst
            title = "Ernst formula"
        elseif typeof(s) == Coherent
            title = string("configuration (n = ", f(:n, s), ")")
        elseif typeof(s) == RFSpoiled
            title = string("RF spoiled, (ϕ = ", f(:ϕ_inc_deg, s), "°)")
        elseif typeof(s) == Balanced
            title = "bSSFP"
        end

        par!(s, Dict(
            :ρ => βs.*ρ,
            :R1 => R1./βs.^2,
            :R2 => R2,
            :TR => TR,
            :β => 1.0,      # set to 1, since we explicitly rely on βs
            :α_nom => αs./βs))
        
        plts[1,is] = plot()

        for (iβ,β) in enumerate(βs)
            plot!(
                rad2deg.(αs[:]),
                reshape(abs.(f(:m, s)[iβ,:]), length(αs)),
                label = string(β),
                margin = margin,
                title = title,
                titlefontsize = 12,
                legend = :outerright, 
                legend_title = "β"
            )
        end
        
        is == length(seqs) && xlabel!(L"$\alpha$ [deg]")
        ylabel!(L"$\beta\cdot |m(R_1/\beta^2,\alpha/\beta)|$")
        
        par!(s, Dict(
            :ρ => ρ,
            :R1 => R1,
            :α_nom => αs))
            
        m_1 = f(:m, s)
        
        par!(s, Dict(
            :ρ => βr.*ρ,
            :R1 => R1./βr.^2,
            :α_nom => αs./βr))
            
        m_β = f(:m, s)

        plts[2,is] = heatmap(rad2deg.(αs[:]),
            βr[:], 
            errfak .* (abs.(m_β./m_1) .- 1), 
            yticks = βs,
            clim = errfak .* errbnds,
            color = ncats == 0 ? cmp : cgrad(cmp, ncats, categorical=true),
            margin = margin,
            title = "rel. deviation [%]",
            titlefontsize = 12)
            
        is == length(seqs) && xlabel!(L"$\alpha$ [deg]")
        ylabel!(L"\beta")
    end

    display(plot(plts..., 
        layout = size(plts'),
        dpi = dpi,
        size = size(plts) .* subplt_size))
end
#= 
==================================================================  
Precision of estimators
==================================================================  
=#

"""
```
∂m_∂θ(seqs, θs)
```
Calculates the matrix `F_ij` of partial derivatives `∂m(s_i)/∂θ_j`
# Parameters
  * `seqs`: Iterable with `typeof(s) ∈ {Ernst, Coherent} for s in seqs`
  * `θs`: Iterable with `typeof(θ) ∈ {:ρ, :β, :R1, :R2} for θ in θs`
"""
∂m_∂θ(seqs, θs) = [∂f(:m, θ, seq) for seq in seqs, θ in θs]

#-----------------------------------------------------------------

"""
```
fisher(seqs, θs, σ)
```
Calculates the Fisher information matrix `∂m_∂θ(seqs, θs)' * ∂m_∂θ(seqs, θs) ./ σ^2`
# Parameters
  * `seqs`: Iterable with `typeof(s) ∈ {Ernst, Coherent} for s in seqs`
  * `θs`: Iterable with `typeof(θ) ∈ {:ρ, :β, :R1, :R2} for θ in θs`
  + `σ`: Noise standard deviation (`length(σ) ∈ {1, length(seqs)}`)
"""
fisher(seqs, θs, σ) = (∂m_∂θ(seqs, θs)' * ∂m_∂θ(seqs, θs)) ./ σ^2

#-----------------------------------------------------------------

"""
```
crb(seqs, θs, σ)
```
Cramer Rao bound (inverse of Fisher matrix): `inv(fisher(seqs, θs, σ))`
# Parameters
  * `seqs`: Iterable with `typeof(s) ∈ {Ernst, Coherent} for s in seqs`
  * `θs`: Iterable with `typeof(θ) ∈ {:ρ, :β, :R1, :R2} for θ in θs`
  + `σ`: Noise standard deviation (`length(σ) ∈ {1, length(seqs)}`)
"""
crb(seqs, θs, σ) = inv(fisher(seqs, θs, σ))

#= 
==================================================================  
SSFP variants: Ernst, Coherent, Balanced
==================================================================  
=#

"""
```
Ernst <: SSFP
```
Implements Ernst formula (cf. supporting information).

Required keys: `(:α_nom, :TR, :ρ, :β, :R1)`
# Constructors
```
Ernst()           # Initiates with empty dictionary.
Ernst(s::SSFP)    # Sequence and tissue parameters are copied from s.
```
"""
mutable struct Ernst <: SSFP
    par::Dict{Symbol,Any}
    Ernst(par) = new(copy(par))
end

Ernst() = Ernst(Dict{Symbol,Any}())
Ernst(s::SSFP) = remove_derived!(Ernst(par(s)))

#-----------------------------------------------------------------

"""
```
Coherent <: Configuration
```
Implements coherent configuration (cf. supporting information).

Required keys: `(:α_nom, :TR, :ρ, :β, :R1, :R2, :n)`
# Constructors
```
Coherent()           # Initiates with empty dictionary.
Coherent(s::SSFP)    # Sequence and tissue parameters are copied from s.
```
"""
mutable struct Coherent <: Configuration             # single configurations
    par::Dict{Symbol,Any}
    Coherent(par) = new(copy(par))
end

Coherent() = Coherent(Dict{Symbol,Any}())
Coherent(s::SSFP) = remove_derived!(Coherent(par(s)))

#-----------------------------------------------------------------

"""
```
RFSpoiled <: Configuration
```
Implements RF spoiled FID 

Higher orders are essentially just more T2* weigthed (n>0)
or negligible (n<0). In case of partial spoiling it is also the 
FID, which is of main interest.

Required keys: `(:α_nom, :TR, :ρ, :β, :R1, :R2, :ϕ_inc)`
# Constructors
```
RFSpoiled()           # Initiates with empty dictionary.
RFSpoiled(s::SSFP)    # Sequence and tissue parameters are copied from s.
```
"""
mutable struct RFSpoiled <: Configuration             # single configurations
    par::Dict{Symbol,Any}
    RFSpoiled(par) = new(copy(par))
end

RFSpoiled() = RFSpoiled(Dict{Symbol,Any}())
RFSpoiled(s::SSFP) = remove_derived!(RFSpoiled(par(s)))

#-----------------------------------------------------------------

"""
```
Balanced <: SSFP
```
Implements balanced SSFP (cf. supporting information).

Required keys: `(:α_nom, :TR, :ρ, :β, :R1, :R2, :ϑ)`
# Constructors
```
Balanced()           # Initiates with empty dictionary.
Balanced(s::SSFP)    # Sequence and tissue parameters are copied from s.
```
"""
mutable struct Balanced <: SSFP             # balanced SSFP
    par::Dict{Symbol,Any}
    Balanced(par) = new(copy(par))
end

Balanced() = Balanced(Dict{Symbol,Any}())
Balanced(s::SSFP) = remove_derived!(Balanced(par(s)))

#= 
==================================================================  
Auxiliary routines (not exported)
==================================================================  
=#

"""
```
α(s::SSFP)
```
Actual flip angle [rad]
"""
α(s::SSFP) = f(:β, s) .* f(:α_nom, s)

#-----------------------------------------------------------------

"""
```
cα(s::SSFP)
```
Cosine of `α(s)`
"""
cα(s::SSFP) = cos.(f(:α, s))

#-----------------------------------------------------------------

"""
```
∂cα(∂x::Symbol, s::SSFP)
```
Partial derivative of `cα(s)` with respect to `∂x ∈ {:β, :α_nom, :b, :cα)` (zero returned otherwise)
"""
function ∂cα(∂x::Symbol, s::SSFP)
    ∂x ∈ (:β, :α_nom) && return - ∂f(:b, ∂x, s)
    ∂x == :b && return - 1.0
    ∂x == :cα ? 1.0 : 0.0
end

#-----------------------------------------------------------------

"""
```
sα(s::SSFP)
```
Sine of `α(s)`
"""
sα(s::SSFP) = sin.(f(:α, s))

#-----------------------------------------------------------------

"""
```
∂sα(∂x::Symbol, s::SSFP)
```
Partial derivative of `sα(s)` with respect to `∂x ∈ {:β, :α_nom, :b, :sα)` (zero returned otherwise)
"""
function ∂sα(∂x::Symbol, s::SSFP)
    ∂x ∈ (:β, :α_nom) && return - ∂f(:sα, :b, s) .* ∂f(:b, ∂x, s)
    ∂x == :b && return (1.0 .- f(:b, s)) ./ f(:sα, s)
    ∂x == :sα ? 1.0 : 0.0
end

#-----------------------------------------------------------------

"""
```
ϕ_inc(s::RFSpoiled)
```
Phase difference increment [rad]
"""
ϕ_inc(s::RFSpoiled) = deg2rad.(round(f(:ϕ_inc_deg, s), digits=f(:rnd_inc, s)))

#-----------------------------------------------------------------

"""
```
y(s::RFSpoiled)
```
Phase factor with respect to RF spoiling: `y = exp(1im * ϕ_inc)`
"""
y(s::RFSpoiled) = exp.(1im .* f(:ϕ_inc, s))

#-----------------------------------------------------------------

"""
```
a(s::SSFP)
```
Function `a`, as defined in the manuscript.
"""
a(s::SSFP) = 1.0 .- f(:E1, s)

#-----------------------------------------------------------------

"""
```
∂a(∂x::Symbol, s::SSFP)
```
Partial derivative of `a(s)` with respect to `∂x ∈ {:R1, :TR, :E1, :a)` (zero returned otherwise)
"""
function ∂a(∂x::Symbol, s::SSFP)
    ∂x ∈ (:R1, :TR) && return - ∂f(:E1, ∂x, s)
    ∂x == :E1 && return - 1.0
    ∂x == :a ? 1.0 : 0.0
end

#-----------------------------------------------------------------

"""
```
b(s::SSFP)
```
Function `b`, as defined in the manuscript.
"""
b(s::SSFP) = 1.0 .- f(:cα, s)

#-----------------------------------------------------------------

"""
```
∂b(∂x::Symbol, s::SSFP)
```
Partial derivative of `b(s)` with respect to `∂x ∈ {:β, :α_nom, :b)` (zero returned otherwise)
"""
function ∂b(∂x::Symbol, s::SSFP)
    ∂x == :β && return f(:α_nom, s) .* sin.(f(:α, s))
    ∂x == :α_nom && return f(:β, s) .* sin.(f(:α, s))
    ∂x == :b ? 1.0 : 0.0
end

#-----------------------------------------------------------------

"""
```
c(s::SSFP)
```
Function `c`, as defined in the manuscript.
"""
c(s::SSFP) = (f(:a, s) .- f(:b, s)) ./ (f(:a, s) .+ f(:b, s) .- f(:a, s) .* f(:b, s))

#-----------------------------------------------------------------

"""
```
∂c(∂x::Symbol, s::SSFP)
```
Partial derivative of `c(s)` with respect to `∂x ∈ {:a, :b, :c)` (zero returned otherwise)
"""
function ∂c(∂x::Symbol, s::SSFP)
    ∂x == :a && return f(:b, s) .* (2.0 .- f(:b, s)) ./
           (f(:a, s) .+ f(:b, s) .- f(:a, s) .* f(:b, s)) .^ 2
    ∂x == :b && return f(:a, s) .* (f(:a, s) .- 2.0) ./
           (f(:a, s) .+ f(:b, s) .- f(:a, s) .* f(:b, s)) .^ 2
    ∂x == :c ? 1.0 : 0.0
end

#-----------------------------------------------------------------

"""
```
d(s::SSFP)
```
Function `d`, as defined in the manuscript.
"""
d(s::SSFP) =
    (
        1.0 .- f(:c, s) .* f(:E2, s) .^ 2 .-
        sqrt.((1.0 .- (f(:c, s) .* f(:E2, s)) .^ 2) .* (1.0 .- f(:E2, s) .^ 2))
    ) ./ (1.0 .- f(:c, s))

#-----------------------------------------------------------------

"""
```
∂d(∂x::Symbol, s::SSFP)
```
Partial derivative of `d(s)` with respect to `∂x ∈ {:a, :b, :c, :E2, :d)` (zero returned otherwise)
"""
function ∂d(∂x::Symbol, s::SSFP)
    ∂x == :a && return ∂f(:d, :c, s) .* ∂f(:c, :a, s)
    ∂x == :b && return ∂f(:d, :c, s) .* ∂f(:c, :b, s)
    sqrt_ = sqrt.((1.0 .- (f(:c, s) .* f(:E2, s)) .^ 2) .* (1.0 .- f(:E2, s) .^ 2))
    if ∂x == :c
        ∂sqrt_ = -f(:c, s) .* f(:E2, s) .^ 2 .* (1.0 .- f(:E2, s) .^ 2) ./ sqrt_
        return (
            -f(:E2, s) .^ 2 - ∂sqrt_ .+
            (1.0 .- f(:c, s) .* f(:E2, s) .^ 2 .- sqrt_) ./ (1.0 .- f(:c, s))
        ) ./ (1.0 .- f(:c, s))
    elseif ∂x == :E2
        ∂sqrt_ =
            -f(:E2, s) .* (1.0 .+ f(:c, s) .^ 2 .* (1.0 .- 2.0 .* f(:E2, s) .^ 2)) ./ sqrt_
        return (-2.0 .* f(:c, s) .* f(:E2, s) - ∂sqrt_) ./ (1.0 .- f(:c, s))
    end
    ∂x == :d ? 1.0 : 0.0
end

#-----------------------------------------------------------------

"""
```
E1(s::SSFP)
```
Longitudinal relaxation: `E1 = exp(- R1 * TR)`
"""
E1(s::SSFP) = exp.(- f(:R1, s) .* f(:TR, s))

#-----------------------------------------------------------------

"""
```
∂E1(∂x::Symbol, s::SSFP)
```
Partial derivative of `E1(s)` with respect to `∂x ∈ {:R1, :TR, :a, :E1)` (zero returned otherwise)
"""
function ∂E1(∂x::Symbol, s::SSFP)
    ∂x == :R1 && return - f(:TR, s) .* f(:E1, s)
    ∂x == :TR && return - f(:R1, s) .* f(:E1, s)
    ∂x == :a && return - 1.0
    ∂x == :E1 ? 1.0 : 0.0
end

#-----------------------------------------------------------------

"""
```
E2(s::SSFP)
```
Transverse relaxation: `E2 = exp(- R2 * TR)`
"""
E2(s::SSFP) = exp.(-f(:R2, s) .* f(:TR, s))

#-----------------------------------------------------------------

"""
```
∂E2(∂x::Symbol, s::SSFP)
```
Partial derivative of `E2(s)` with respect to `∂x ∈ {:R2, :TR, :E2)` (zero returned otherwise)
"""
function ∂E2(∂x::Symbol, s::SSFP)
    ∂x == :R2 && return -f(:TR, s) .* f(:E2, s)
    ∂x == :TR && return -f(:R2, s) .* f(:E2, s)
    ∂x == :E2 ? 1.0 : 0.0
end

#-----------------------------------------------------------------

"""
```
A(s::SSFP)
```
Function `A`, as defined in the manuscript.
"""
A(s::SSFP) = f(:E2, s) .* (1 .+ f(:c, s)) ./ (2.0 .- (1.0 .- f(:c, s)) .* f(:d, s))

#-----------------------------------------------------------------

"""
```
∂A(∂x::Symbol, s::SSFP)
```
Partial derivative of `A(s)` with respect to `∂x ∈ {:a, :b, :c, :d, :E2, :A)` (zero returned otherwise)
"""
function ∂A(∂x::Symbol, s::SSFP)
    ∂x ∈ (:a, :b) && return ∂f(:A, :c, s) .* ∂f(:c, ∂x, s) .+ ∂f(:A, :d, s) .* ∂f(:d, ∂x, s)
    den = 2.0 .- (1.0 .- f(:c, s)) .* f(:d, s)
    ∂x == :c && return f(:E2, s) .* (1.0 .- (1.0 .+ f(:c, s)) .* f(:d, s) ./ den) ./ den
    ∂x == :d && return f(:E2, s) .* (1.0 .+ f(:c, s)) .* (1.0 .- f(:c, s)) ./ (den .* den)
    ∂x == :E2 && return (1 .+ f(:c, s)) .*
           (1.0 .+ f(:E2, s) .* (1.0 .- f(:c, s)) .* ∂f(:d, :E2, s) ./ den) ./ den
    ∂x == :A ? 1.0 : 0.0
end

#-----------------------------------------------------------------

"""
```
m(s::Ernst)
```
Ernst formula (= steady state for ideal spoiling)
"""
m(s::Ernst) =
    f(:ρ, s) .* f(:sα, s) .* f(:a, s) ./ (f(:a, s) .+ f(:b, s) .- f(:a, s) .* f(:b, s))

#-----------------------------------------------------------------

"""
```
∂m(∂x::Symbol, s::Ernst)
```
Partial derivative of Ernst formula with respect to `∂x ∈ {:a, :b, :E2}` (throws ArgumentError otherwise)
"""
function ∂m(∂x::Symbol, s::Ernst)
    ∂x ∉ (:a, :b, :E2) && throw(ArgumentError("unrecognized derivative"))
    den = f(:a, s) .+ f(:b, s) .- f(:a, s) .* f(:b, s)
    ∂x == :a && return f(:ρ, s) .* f(:sα, s) .* f(:b, s) ./ den .^ 2
    ∂x == :b && return f(:ρ, s) .*
           (∂f(:sα, :b, s) .+ f(:sα, s) .* (f(:a, s) .- 1.0) ./ den) .*
           f(:a, s) ./ den
    ∂x == :E2 && return 0.0
end

#-----------------------------------------------------------------

"""
```
m(s::Coherent)
```
Steady-State of coherent configuration
"""
function m(s::Coherent)
    er = Ernst(s)
    m_ = f(:m, er) ./ (1.0 .+ f(:c, s) .* f(:d, s))
    f(:n, s) < 0 && (m_ = -m_ .* f(:d, s) ./ f(:E2, s))
    f(:n, s) ∉ (0, -1) &&
        (m_ = m_ .* f(:A, s) .^ (f(:n, s) >= 0 ? f(:n, s) : -(1 + f(:n, s))))
    m_
end

#-----------------------------------------------------------------

"""
```
∂m(∂x::Symbol, s::Coherent)
```
Partial derivative of coherent SSFP configuration with respect to `∂x ∈ {:a, :b, :E2}` (throws ArgumentError otherwise)
"""
function ∂m(∂x::Symbol, s::Coherent)
    ∂x ∉ (:a, :b, :E2) && throw(ArgumentError("unrecognized derivative"))

    er = Ernst(s)

    (m_, ∂m_) = (
        f(:m, er) ./ (1.0 .+ f(:c, s) .* f(:d, s)),
        (
            ∂f(:m, ∂x, er) .-
            f(:m, er) .* (∂f(:c, ∂x, s) .* f(:d, s) .+ f(:c, s) .* ∂f(:d, ∂x, s)) ./
            (1.0 .+ f(:c, s) .* f(:d, s))
        ) ./ (1.0 .+ f(:c, s) .* f(:d, s)),
    )

    f(:n, s) < 0 && (
        (m_, ∂m_) = (
            -m_ .* f(:d, s) ./ f(:E2, s),
            -∂m_ .* f(:d, s) ./ f(:E2, s) .-
            m_ .* (∂f(:d, ∂x, s) .- f(:d, s) .* ∂f(:E2, ∂x, s) ./ f(:E2, s)) ./ f(:E2, s),
        )
    )

    (n_exp = f(:n, s) >= 0 ? f(:n, s) : -(1 + f(:n, s))) != 0 && (
        ∂m_ =
            ∂m_ .* f(:A, s) .^ n_exp .+
            n_exp .* m_ .* f(:A, s) .^ (n_exp - 1) .* ∂f(:A, ∂x, s)
    )

    return ∂m_
end

#-----------------------------------------------------------------

"""
```
m(s::Balanced)
```
Steady-state of balanced SSFP
"""
function m(s::Balanced)
    fid, echo, eiϑ =
        par!(Coherent(s), :n, 0), par!(Coherent(s), :n, -1), exp.(im .* f(:ϑ, s))
    f(:m, fid) ./ (1.0 .- f(:A, fid) .* eiϑ) .+ f(:m, echo) ./ (eiϑ .- f(:A, fid))
end

#-----------------------------------------------------------------

"""
```
∂m(∂x::Symbol, s::Balanced)
```
Partial derivative of balanced SSFP with respect to `∂x ∈ {:a, :b, :E2}` (throws ArgumentError otherwise)
"""
function ∂m(∂x::Symbol, s::Balanced)
    ∂x ∉ (:a, :b, :E2) && throw(ArgumentError("unrecognized derivative"))
    fid, echo, eiϑ =
        par!(Coherent(s), :n, 0), par!(Coherent(s), :n, -1), exp.(im .* f(:ϑ, s))
    (
        ∂f(:m, ∂x, fid) .+
        f(:m, fid) .* ∂f(:A, ∂x, fid) .* eiϑ ./ (1.0 .- f(:A, fid) .* eiϑ)
    ) ./ (1.0 .- f(:A, fid) .* eiϑ) .+
    (∂f(:m, ∂x, echo) .+ f(:m, echo) .* ∂f(:A, ∂x, fid) ./ (eiϑ .- f(:A, fid))) ./
    (eiϑ .- f(:A, fid))
end

#-----------------------------------------------------------------

"""
```
rel_err(::RFSpoiled)
```
Default relative accuracy, if not set explicitly
"""
rel_err(::RFSpoiled) = eps()

#-----------------------------------------------------------------

"""
```
rnd_inc(::RFSpoiled)
```
Number of supported digits in the specification of `ϕ_inc_deg`.
Required to determine (and limit) the periodicity.
"""
rnd_inc(::RFSpoiled) = 1

#-----------------------------------------------------------------

"""
```
max_iter(::RFSpoiled)
```
Number of allowed iterations.
"""
max_iter(::RFSpoiled) = 1000

#-----------------------------------------------------------------

"""
```
Λ(s::RFSpoiled)
```
Λ, as defined by the continuous fraction expansion (27) in Ganter, MRM (2006) 55:98-107
"""
function Λ(s::RFSpoiled)
    # to avoid conflicts with occupied dimensions,
    # we need some reshaping such that powers of y = exp(iϕ) (see below)
    # can be stored in an extra dimension

    (cα_, sα_, E1_, E2_) = (f(:cα,s), f(:sα,s), f(:E1,s), f(:E2,s))
    z_ = zeros(size(cα_)) .+ zeros(size(E1_)) .+ zeros(size(E2_))
    if size(z_) != ()
        cα_ = (cα_ .+ z_)[:]
        sα_ = (sα_ .+ z_)[:]
        E1_ = (E1_ .+ z_)[:]
        E2_ = (E2_ .+ z_)[:]
    end
    
    # avoid divisions by zero
    
    tiny_num = eps()    

    # since we limit the precision of ϕ_inc_deg, we can exploit periodicity (defined by N)

    fak = 10^round(Int, f(:rnd_inc, s))
    N = 360fak ÷ gcd(360fak, round(Int, fak * f(:ϕ_inc_deg, s)))
    N == 1 && (N = 2)
    
    ns = reshape(collect(1:N), 1, N)
    yn = f(:y,s) .^ ns
    yn2 = yn .* yn
    
    αn = 1 .- cα_ .* E1_ .* yn
    βn = cα_ .- E1_ .* yn
    γn = 0.5(αn .- βn)
    
    γα = γn ./ αn
    βγ = βn ./ γn
   
    E22y = E2_.^2 .* f(:y,s)
    
    a0 = (γα[:,1] .+ βγ[:,1]) .* E22y
    b0 = - βγ[:,1] .* E22y
   
    n1 = mod.(ns, N)[:] .+ 1
   
    an = - γα .* (γα[:,n1] .+ βγ[:,n1]) .* yn2 .* E22y
    bn = 1 .+ γα .* βγ[:,n1] .* yn2 .* E22y
    
    Λ_ = b0
    
    C0 = Λ_
    D1 = bn[:,1]
    D1[D1 .== 0] .= tiny_num
    
    C1 = bn[:,1] .+ a0 ./ C0
    C1[C1 .== 0] .= tiny_num
    
    D1 = 1 ./ D1
    dj = C1 .* D1
    
    Λ_ = Λ_ .* dj
    
    C0 = C1
    D0 = D1
    
    j = 1
    j1 = 2
    
    iter = 0
    
    while max(abs.(dj .- 1)...) > f(:rel_err,s) && iter < f(:max_iter,s)

        D1 = bn[:,j1] .+ an[:,j] .* D0
        D1[D1 .== 0] .= tiny_num
        
        C1 = bn[:,j1] .+ an[:,j] ./ C0
        C1[C1 .== 0] .= tiny_num
        
        D1 = 1 ./ D1
        dj = C1 .* D1
        
        Λ_ = Λ_ .* dj
        
        C0 = C1
        D0 = D1
        
        j = mod(j, N) + 1
        j1 = mod(j, N) + 1
        
        iter = iter + 1        
        
    end
    
    return size(z_) == () ? Λ_[] : reshape(Λ_, size(z_))
end

#-----------------------------------------------------------------

"""
```
∂Λ(∂x::Symbol, s::RFSpoiled)
```
Partial derivative of `Λ(s::RFSpoiled)` with respect to `∂x ∈ {:a, :b, :E2)` (throws ArgumentError otherwise)
"""
function ∂Λ(∂x, s::RFSpoiled)
    ∂x ∉ (:a, :b, :E2) && throw(ArgumentError("unrecognized derivative"))
    
    # to avoid conflicts with occupied dimensions,
    # we need some reshaping such that powers of y = exp(iϕ) (see below)
    # can be stored in an extra dimension

    (cα_, sα_, E1_, E2_) = (f(:cα,s), f(:sα,s), f(:E1,s), f(:E2,s))
    (∂cα_, ∂sα_, ∂E1_, ∂E2_) = (∂f(:cα, ∂x, s), ∂f(:sα, ∂x, s), ∂f(:E1, ∂x, s), ∂f(:E2, ∂x, s))
    z_ = zeros(size(cα_)) .+ zeros(size(E1_)) .+ zeros(size(E2_))
    if size(z_) != ()
        cα_ = (cα_ .+ z_)[:]
        sα_ = (sα_ .+ z_)[:]
        E1_ = (E1_ .+ z_)[:]
        E2_ = (E2_ .+ z_)[:]
        ∂cα_ = (∂cα_ .+ z_)[:]
        ∂sα_ = (∂sα_ .+ z_)[:]
        ∂E1_ = (∂E1_ .+ z_)[:]
        ∂E2_ = (∂E2_ .+ z_)[:]
    end
    
    # avoid divisions by zero
    
    tiny_num = eps()    

    # since we limit the precision of ϕ_inc_deg, we can exploit periodicity (defined by N)

    fak = 10^round(Int, f(:rnd_inc, s))
    N = 360fak ÷ gcd(360fak, round(Int, fak * f(:ϕ_inc_deg, s)))
    N == 1 && (N = 2)

    ns = reshape(collect(1:N), 1, N)
    yn = f(:y,s) .^ ns
    yn2 = yn .* yn
    
    αn = 1 .- cα_ .* E1_ .* yn
    βn = cα_ .- E1_ .* yn
    γn = 0.5(αn .- βn)
    ∂αn = - (∂cα_ .* E1_ .+ cα_ .* ∂E1_) .* yn
    ∂βn = ∂cα_ .- ∂E1_ .* yn
    ∂γn = 0.5(∂αn .- ∂βn)

    γα = γn ./ αn
    βγ = βn ./ γn
    ∂γα = (∂γn .- ∂αn .* γn ./ αn) ./ αn
    ∂βγ = (∂βn .- ∂γn .* βn ./ γn) ./ γn
    
    E22y = E2_.^2 .* f(:y,s)
    ∂E22y = 2 .* E2_ .* ∂E2_ .* f(:y,s)
    
    a0 = (γα[:,1] .+ βγ[:,1]) .* E22y
    b0 = - βγ[:,1] .* E22y
    ∂a0 = (∂γα[:,1] .+ ∂βγ[:,1]) .* E22y .+ (γα[:,1] .+ βγ[:,1]) .* ∂E22y
    ∂b0 = - ∂βγ[:,1] .* E22y .- βγ[:,1] .* ∂E22y
    
    n1 = mod.(ns, N)[:] .+ 1
   
    an = - γα .* (γα[:,n1] .+ βγ[:,n1]) .* yn2 .* E22y
    bn = 1 .+ γα .* βγ[:,n1] .* yn2 .* E22y
    ∂an = - ∂γα .* (γα[:,n1] .+ βγ[:,n1]) .* yn2 .* E22y .-
        γα .* (∂γα[:,n1] .+ ∂βγ[:,n1]) .* yn2 .* E22y .-
        γα .* (γα[:,n1] .+ βγ[:,n1]) .* yn2 .* ∂E22y
    ∂bn = ∂γα .* βγ[:,n1] .* yn2 .* E22y .+
        γα .* ∂βγ[:,n1] .* yn2 .* E22y .+
        γα .* βγ[:,n1] .* yn2 .* ∂E22y

    Λ_ = b0
    ∂Λ_ = ∂b0
    
    C0 = Λ_
    D1 = bn[:,1]
    D1[D1 .== 0] .= tiny_num
    ∂C0 = ∂Λ_
    ∂D1 = ∂bn[:,1]
    
    C1 = bn[:,1] .+ a0 ./ C0
    C1[C1 .== 0] .= tiny_num
    ∂C1 = ∂bn[:,1] .+ (∂a0 .- ∂C0 .* a0 ./ C0) ./ C0
    
    ∂D1 = - ∂D1 ./ D1.^2
    D1 = 1 ./ D1
    dj = C1 .* D1
    ∂dj = ∂C1 .* D1 .+ C1 .* ∂D1
    
    ∂Λ_ = ∂Λ_ .* dj .+ Λ_ .* ∂dj
    Λ_ = Λ_ .* dj
    
    C0 = C1
    D0 = D1
    ∂C0 = ∂C1
    ∂D0 = ∂D1
    
    j = 1
    j1 = 2
    
    iter = 0
    
    while max(abs.(dj .- 1)...) > f(:rel_err,s) && iter < f(:max_iter,s)

        D1 = bn[:,j1] .+ an[:,j] .* D0
        D1[D1 .== 0] .= tiny_num
        ∂D1 = ∂bn[:,j1] .+ ∂an[:,j] .* D0 .+ an[:,j] .* ∂D0 
        
        C1 = bn[:,j1] .+ an[:,j] ./ C0
        C1[C1 .== 0] .= tiny_num
        ∂C1 = ∂bn[:,j1] .+ (∂an[:,j] .- ∂C0 .* an[:,j] ./ C0) ./ C0
        
        ∂D1 = - ∂D1 ./ D1.^2
        D1 = 1 ./ D1
        dj = C1 .* D1
        ∂dj = ∂C1 .* D1 .+ C1 .* ∂D1
        
        ∂Λ_ = ∂Λ_ .* dj .+ Λ_ .* ∂dj
        Λ_ = Λ_ .* dj
        
        C0 = C1
        D0 = D1
        ∂C0 = ∂C1
        ∂D0 = ∂D1
        
        j = mod(j, N) + 1
        j1 = mod(j, N) + 1
        
        iter = iter + 1        
        
    end
    
    return size(z_) == () ? ∂Λ_[] : reshape(∂Λ_, size(z_))
end

#-----------------------------------------------------------------

"""
```
nom(s::RFSpoiled)
```
Nominator of RF spoiled SSFP
"""
nom(s::RFSpoiled) = f(:sα,s) .* (1 .- conj.(f(:Λ,s))) .* f(:a,s)

#-----------------------------------------------------------------

"""
```
∂nom(∂x::Symbol, s::RFSpoiled)
```
Partial derivative of `nom(s::RFSpoiled)` with respect to `∂x ∈ {:a, :b, :E2)` (throws ArgumentError otherwise)
"""
function ∂nom(∂x::Symbol, s::RFSpoiled)
    ∂x ∉ (:a, :b, :E2) && throw(ArgumentError("unrecognized derivative"))
    
    return ∂f(:sα, ∂x, s) .* (1 .- conj.(f(:Λ,s))) .* f(:a,s) .-
        f(:sα,s) .* conj.(∂f(:Λ, ∂x, s)) .* f(:a,s) .+
        f(:sα,s) .* (1 .- conj.(f(:Λ,s))) .* ∂f(:a, ∂x, s)
end

#-----------------------------------------------------------------

"""
```
den(s::RFSpoiled)
```
Denominator of RF spoiled SSFP
"""
den(s::RFSpoiled) = f(:a, s) .+ 
    f(:b, s) .* (1 .- f(:a, s) .- (2 .- f(:a,s)) .* real.(f(:Λ,s))) .+
    (f(:b,s) .- f(:a,s)) .* abs.(f(:Λ,s)).^2

#-----------------------------------------------------------------

"""
```
∂den(∂x::Symbol, s::RFSpoiled)
```
Partial derivative of `den(s::RFSpoiled)` with respect to `∂x ∈ {:a, :b, :E2)` (throws ArgumentError otherwise)
"""
function ∂den(∂x::Symbol, s::RFSpoiled)
    ∂x ∉ (:a, :b, :E2) && throw(ArgumentError("unrecognized derivative"))
    
    return ∂f(:a, ∂x, s) .+ 
    ∂f(:b, ∂x, s) .* (1 .- f(:a, s) .- (2 .- f(:a,s)) .* real.(f(:Λ,s))) .+
    f(:b, s) .* (- ∂f(:a, ∂x, s) .+ ∂f(:a, ∂x, s) .* real.(f(:Λ,s)) .- (2 .- f(:a,s)) .* real.(∂f(:Λ, ∂x, s))) .+
    (∂f(:b, ∂x, s) .- ∂f(:a, ∂x, s)) .* abs.(f(:Λ,s)).^2 .+
    (f(:b,s) .- f(:a,s)) .* 2 .* real.(∂f(:Λ, ∂x, s) .* conj.(f(:Λ,s)))
end

#-----------------------------------------------------------------

"""
```
m(s::RFSpoiled)
```
Steady-state of RF spoiled SSFP (magnitude only!)
"""
m(s::RFSpoiled) = abs.(f(:ρ,s) .* f(:nom,s) ./ f(:den,s))

#-----------------------------------------------------------------

"""
```
∂m(::Symbol, ::RFSpoiled)
```
Partial derivative of RF spoiled SSFP with respect to `∂x ∈ {:a, :b, :E2}` (throws ArgumentError otherwise)
"""
function ∂m(∂x::Symbol, s::RFSpoiled)
    ∂x ∉ (:a, :b, :E2) && throw(ArgumentError("unrecognized derivative"))
    
    z = f(:ρ,s) .* f(:nom,s) ./ f(:den,s)
    ∂z = f(:ρ,s) .* (∂f(:nom, ∂x, s) .- ∂f(:den, ∂x, s) .* f(:nom,s) ./ f(:den,s)) ./ f(:den,s)
    return real.(∂z .* conj.(z)) ./ f(:m,s)
end

end  # module
