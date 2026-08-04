# Does the small-γ mass-derivative cancellation survive to a whole trajectory,
# at timesteps anyone would actually use?
#
# Run:
#   julia --project=perf perf/nbody-gradient-endtoend.jl
#
# `nbody-mass-derivative.jl` showed the ∂/∂m_i column of one pair kernel loses
# ~4.5 digits at γ = 0.01. That is a statement about one sub-step in isolation.
# What matters is the gradient of a whole integration, and whether γ = 0.01 is
# a place anyone stands: γ per step ≈ √(GM/a)·h/r, so the recommended
# h ≈ P_min/20 puts a circular pair at γ ≈ 0.31, and γ = 0.01 needs h ≈ P/600.
#
# This sweeps h across that range on a Kepler-36-like three-body system and
# reports, against BigFloat(256) evaluation of identical expressions:
#   - d(final state)/dm  — the mass gradient the cancellation threatens
#   - d(final state)/dh  — the quantity upstream's `dTime` path is flagged as
#     getting wrong ("[Currently this routine is not giving the correct dqdt
#     values. -EA 8/12/2019]"), which we get from the same Duals for free.
# Errors are relative to the gradient norm, never element-by-element.

using StaticArrays
using ForwardDiff
using Printf
using PlanetOrbits: NBodyState, ahl21_step, G_au3_day2, _delxv_gamma, _comp_sum, _dot3

setprecision(BigFloat, 256)

include(joinpath(@__DIR__, "nbody-ref-step.jl"))

# Kepler-36-like: two close-in super-Earths near a 7:6 period ratio.
const MASSES = (1.071, 1.36e-5, 2.34e-5)     # M⊙
const SMA    = (0.0, 0.1153, 0.1283)         # AU
const PHASE  = (0.0, 0.4, 2.1)

function initial_state(::Type{T}, m) where {T}
    N = 3
    x = zeros(T, 3, N); v = zeros(T, 3, N)
    for i in 2:N
        a = T(SMA[i])
        vc = sqrt(T(G_au3_day2) * m[1] / a)
        ph = T(PHASE[i])
        x[1, i] = a * cos(ph); x[2, i] = a * sin(ph); x[3, i] = a * T(0.002)
        v[1, i] = -vc * sin(ph); v[2, i] = vc * cos(ph)
    end
    for k in 1:3                              # zero the barycentre
        x[k, 1] = -sum(m[i] * x[k, i] for i in 2:N) / m[1]
        v[k, 1] = -sum(m[i] * v[k, i] for i in 2:N) / m[1]
    end
    return NBodyState(SMatrix{3,N,T}(x), SMatrix{3,N,T}(v))
end

# Integrate nstep steps of size h and return the flattened final state.
function integrate(m::SVector{3}, h, nstep::Int)
    T = promote_type(eltype(m), typeof(h))
    Gm = SVector{3,T}(T(G_au3_day2) .* m)
    st = initial_state(T, SVector{3,T}(m))
    hT = T(h)
    for _ in 1:nstep
        st = ahl21_step(st, Gm, hT)
    end
    return SVector{18,T}(ntuple(l -> begin
        k = (l - 1) % 3 + 1; j = (l - 1) ÷ 3 + 1
        j <= 3 ? st.x[k, j] + st.xerr[k, j] : st.v[k, j - 3] + st.verr[k, j - 3]
    end, Val(18)))
end

# Same trajectory through the Matrix-based reference map, so BigFloat and
# Dual{BigFloat} can run it.
function integrate_ref(m::AbstractVector, h, nstep::Int)
    T = promote_type(eltype(m), typeof(h))
    Gm = T.(T(G_au3_day2) .* m)
    st = initial_state(T, SVector{3,T}(m))
    x = Matrix{T}(st.x); v = Matrix{T}(st.v)
    xerr = zeros(T, 3, 3); verr = zeros(T, 3, 3)
    for _ in 1:nstep
        ref_step(x, v, xerr, verr, collect(Gm), T(h))
    end
    return SVector{18,T}(ntuple(l -> begin
        k = (l - 1) % 3 + 1; j = (l - 1) ÷ 3 + 1
        j <= 3 ? x[k, j] + xerr[k, j] : v[k, j - 3] + verr[k, j - 3]
    end, Val(18)))
end

# γ of the inner pair over one step, for orientation.
inner_gamma(h) = h * sqrt(G_au3_day2 * MASSES[1] / SMA[2]) / SMA[2]

const P_INNER = 2π * sqrt(SMA[2]^3 / (G_au3_day2 * MASSES[1]))

function run()
    @printf("Kepler-36-like 3-body. Inner period %.3f d. Errors vs BigFloat(256),\n", P_INNER)
    println("relative to the gradient norm. `∂/∂m` is the 3×18 mass Jacobian,")
    println("`∂/∂h` is the step-size derivative (upstream's flagged dqdt).\n")
    @printf("%10s %8s %8s %6s | %11s %11s | %9s\n",
            "h", "h/P_inner", "γ/step", "nstep", "∂/∂m relerr", "∂/∂h relerr", "ref==ahl21")
    println("-"^84)
    for div in (20, 60, 200, 600, 2000)
        h = P_INNER / div
        nstep = max(4, round(Int, 3 * P_INNER / h))     # ~3 inner orbits either way
        m64 = SVector{3,Float64}(MASSES)
        mbig = SVector{3,BigFloat}(big.(MASSES))

        # Both sides through the same reference map, so this isolates precision.
        Jm64 = ForwardDiff.jacobian(m -> integrate_ref(m, h, nstep), m64)
        Jmbig = ForwardDiff.jacobian(m -> integrate_ref(m, big(h), nstep), mbig)
        dh64 = ForwardDiff.derivative(hh -> integrate_ref(m64, hh, nstep), h)
        dhbig = ForwardDiff.derivative(hh -> integrate_ref(mbig, hh, nstep), big(h))

        em = Float64(maximum(abs, BigFloat.(Jm64) .- Jmbig) / maximum(abs, Jmbig))
        eh = Float64(maximum(abs, BigFloat.(dh64) .- dhbig) / maximum(abs, dhbig))
        # The reference map must reproduce the shipped kernel in Float64.
        agree = maximum(abs, integrate(m64, h, nstep) .- integrate_ref(m64, h, nstep)) /
                maximum(abs, integrate(m64, h, nstep))
        @printf("%10.5f %8.4f %8.4f %6d | %11.2e %11.2e | %9.1e\n",
                h, h / P_INNER, inner_gamma(h), nstep, em, eh, agree)
    end
end

run()
