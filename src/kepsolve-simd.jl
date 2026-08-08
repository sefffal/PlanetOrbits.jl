# ---------------------------------------------------
# SIMD-batched Kepler solving (Float64 value path)
#
# Vectorizable kernels for solving Kepler's equation across a whole epoch
# column at once: no libm calls, no throw branches — everything reduces to
# FP arithmetic and selects so the LLVM loop vectorizer can pack lanes
# (4-wide on AVX2). Validated against the scalar Markley solver to ≤4e-15
# over M ∈ [-40, 40], e ∈ [0, 0.95].
#
# Bound orbits only. ForwardDiff Duals reach the same kernels through
# `solve_row_dual_simd!` below, which solves the primal roots here and attaches
# partials with the implicit rule from kepsolve.jl.
# ---------------------------------------------------

const PI2_HI = 6.283185307179586
const PI2_LO = 2.4492935982947064e-16
const PIO2_HI = 1.5707963267948966
const PIO2_LO = 6.123233995736766e-17

# rem2pi(x, RoundNearest), branch-free, fma-compensated product.
@inline function vrem2pi(x::Float64)
    rn = round(x * (1 / PI2_HI))
    p = rn * PI2_HI
    perr = fma(rn, PI2_HI, -p)
    return ((x - p) - perr) - rn * PI2_LO
end

# sincos via quadrant reduction + polynomial kernels (|r| ≤ π/4 + small;
# truncation error below ~5e-17 at π/4). Branch-free: selects only.
@inline function _sin_kernel(r::Float64)
    z = r * r
    p = @evalpoly(z, 1.0, -0.16666666666666666, 0.008333333333333333,
        -0.0001984126984126984, 2.7557319223985893e-6, -2.505210838544172e-8,
        1.6059043836821613e-10, -7.647163731819816e-13)
    return r * p
end
@inline function _cos_kernel(r::Float64)
    z = r * r
    return @evalpoly(z, 1.0, -0.5, 0.041666666666666664, -0.001388888888888889,
        2.48015873015873e-5, -2.755731922398589e-7, 2.08767569878681e-9,
        -1.1470745597729725e-11)
end
@inline function vsincos(x::Float64)
    k = round(x * (2 / π))
    r = (x - k * PIO2_HI) - k * PIO2_LO
    q = unsafe_trunc(Int64, k) & 3   # quadrant
    s = _sin_kernel(r)
    c = _cos_kernel(r)
    sinx = ifelse(q == 0, s, ifelse(q == 1, c, ifelse(q == 2, -s, -c)))
    cosx = ifelse(q == 0, c, ifelse(q == 1, -s, ifelse(q == 2, -c, s)))
    return sinx, cosx
end

# sincos for the tiny Markley correction δ (|δ| ≪ 1): Taylor, 3 terms.
@inline function vsincos_small(δ::Float64)
    z = δ * δ
    s = δ * @evalpoly(z, 1.0, -0.16666666666666666, 0.008333333333333333)
    c = @evalpoly(z, 1.0, -0.5, 0.041666666666666664)
    return s, c
end

# Branch-free cbrt for x ≥ 0: exponent-hack seed + 4 Newton steps
# (~1e-14 relative; ample for the Markley E1 seed, whose δ3..δ5 corrections
# converge to the true root at 5th order regardless).
@inline function vcbrt(x::Float64)
    u = reinterpret(UInt64, x)
    y = reinterpret(Float64, u ÷ 3 + 0x2a9f5cc62cb0f9e1)
    y = (2 * y + x / (y * y)) * (1 / 3)
    y = (2 * y + x / (y * y)) * (1 / 3)
    y = (2 * y + x / (y * y)) * (1 / 3)
    y = (2 * y + x / (y * y)) * (1 / 3)
    return ifelse(x == 0, 0.0, y)
end

# Markley (1995) without the early-return branches, on vectorizable kernels.
# Returns sincos(E) directly — the only quantities the observables need — so
# exactly ONE full sincos is evaluated per solve.
@inline function markley_sincosE(Min::Float64, e::Float64)
    M = vrem2pi(Min)
    pi2 = π * π
    α = (3 * pi2 + 8 * (pi2 - π * abs(M)) / (5 * (1 + e))) / (pi2 - 6)
    d = 3 * (1 - e) + α * e
    q = 2 * α * d * (1 - e) - M * M
    r = 3 * α * d * (d - 1 + e) * M + M * M * M
    # abs() inside sqrt: the discriminant is provably ≥ 0 for 0 ≤ e < 1, but
    # without it LLVM keeps a domain-error throw branch that blocks the
    # loop vectorizer.
    w = vcbrt(abs2(abs(r) + sqrt(abs(q * q * q + r * r))))
    E1 = (2 * r * w / @evalpoly(w, q * q, q, 1) + M) / d
    sE1, cE1 = vsincos(E1)
    f2 = e * sE1
    f3 = e * cE1
    f0 = E1 - f2 - M
    f1 = 1 - f3
    δ3 = -f0 / (f1 - f0 * f2 / (2 * f1))
    δ4 = -f0 / @evalpoly(δ3, f1, f2 / 2, f3 / 6)
    δ5 = -f0 / @evalpoly(δ4, f1, f2 / 2, f3 / 6, -f2 / 24)
    # sincos(E1 + δ5) via small-angle rotation — no second full sincos
    sδ, cδ = vsincos_small(δ5)
    sE = sE1 * cδ + cE1 * sδ
    cE = cE1 * cδ - sE1 * sδ
    return sE, cE
end

# SIMD Float64 batch path for one hierarchy row.
function solve_row_simd!(traj::Trajectory{Float64}, row::Row{Float64}, j::Int)
    n_per_day = row.n / year2day_julian
    t_em = traj.t_em
    rx = traj.rx; ry = traj.ry; rz = traj.rz
    rvx = traj.rvx; rvy = traj.rvy; rvz = traj.rvz
    e = row.e
    @inbounds @simd ivdep for k in eachindex(t_em)
        MA = n_per_day * (t_em[k] - row.tp)
        sE, cE = markley_sincosE(MA, e)
        x, y, z, vx, vy, vz = _states_from_E(row, sE, cE)
        rx[k, j] = x; ry[k, j] = y; rz[k, j] = z
        rvx[k, j] = vx; rvy[k, j] = vy; rvz[k, j] = vz
    end
    return traj
end

# ---------------------------------------------------
# The same batch, under ForwardDiff Duals.
#
# The implicit rule in kepsolve.jl already confines the Kepler solve itself to
# primal values — but it reached them through the *scalar* branchy Markley
# solver, one epoch at a time, so a Dual-valued row lost the epoch batching the
# Float64 path gets. That is what design §10.4.1 measured as `propagate!`
# jumping 4.07 -> 23.84 µs at the very *first* partial, long before partial
# arithmetic can account for it.
#
# Nothing in the rule requires the scalar solver. The primal solve is pure
# Float64 work over a contiguous epoch range, so it runs through the same
# `markley_sincosE` kernel the value path uses, and `_dual_sincosE` attaches
# the partials afterwards. Two consequences beyond speed: no transcendental is
# evaluated in Dual arithmetic at all, and the value and gradient paths now
# solve Kepler's equation with the *identical* kernel, so the value carried
# through a gradient evaluation is bit-identical to a plain Float64 evaluation
# rather than agreeing to ~4e-15.
#
# The batching is a fully-unrolled block rather than a `@simd` loop because the
# primal values sit at stride N+1 inside the Dual column, which the loop
# vectorizer will not gather across. A straight-line run of independent,
# branch-free solves is what LLVM's SLP pass packs instead.
# ---------------------------------------------------

# Epochs per straight-line block. 8 keeps two AVX2 vectors in flight through
# the solve without spilling the ~20 live values `markley_sincosE` carries;
# see `perf/dual-passes.jl`, which sweeps it.
const DUAL_SIMD_BLOCK = 8

@inline function _store_row_dual!(traj::Trajectory, row::Row, j::Int, k::Int,
                                  MA::Dual, sE, cE, e::Dual)
    sED, cED = _dual_sincosE(MA, e, sE, cE)
    x, y, z, vx, vy, vz = _states_from_E(row, sED, cED)
    @inbounds begin
        traj.rx[k, j] = x; traj.ry[k, j] = y; traj.rz[k, j] = z
        traj.rvx[k, j] = vx; traj.rvy[k, j] = vy; traj.rvz[k, j] = vz
    end
    return nothing
end

# Unroll `f(1) … f(I)`. Each step is its own method instantiation, so the block
# is straight-line code by construction rather than by the loop unroller's
# discretion — which is the whole point here.
@inline _unroll(f, ::Val{0}) = nothing
@inline function _unroll(f, ::Val{I}) where {I}
    _unroll(f, Val(I - 1))
    f(I)
    return nothing
end

solve_row_dual_simd!(traj::Trajectory, row::Row, j::Int) =
    solve_row_dual_simd!(traj, row, j, Val(DUAL_SIMD_BLOCK))

function solve_row_dual_simd!(traj::Trajectory{D}, row::Row, j::Int,
                              ::Val{B}) where {Tg,N,D<:Dual{Tg,Float64,N},B}
    n_per_day = row.n / year2day_julian
    tp = row.tp
    t_em = traj.t_em
    # A Float64 row lifts to zero partials, which keeps one kernel for the
    # mixed case (plain elements, differentiated frame) at the cost of N
    # provably-zero multiply-adds per epoch.
    e = convert(D, row.e)
    ev = value(e)
    nk = length(t_em)
    k = 1
    @inbounds while k + B - 1 <= nk
        k0 = k
        MA = ntuple(i -> n_per_day * (t_em[k0+i-1] - tp), Val(B))
        # One straight-line block of independent, branch-free primal solves.
        sc = ntuple(i -> markley_sincosE(value(MA[i]), ev), Val(B))
        _unroll(i -> _store_row_dual!(traj, row, j, k0 + i - 1,
                                      MA[i], sc[i][1], sc[i][2], e), Val(B))
        k += B
    end
    @inbounds while k <= nk
        MA = n_per_day * (t_em[k] - tp)
        sE, cE = markley_sincosE(value(MA), ev)
        _store_row_dual!(traj, row, j, k, MA, sE, cE, e)
        k += 1
    end
    return traj
end
