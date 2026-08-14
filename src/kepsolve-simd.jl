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

"""
    PlanetOrbits.VREM2PI_MAX

Largest `|x|` for which `vrem2pi` is a valid stand-in for
`Base.rem2pi(x, RoundNearest)`. **2.0^48 ≈ 2.8e14.**

`vrem2pi` is one Cody–Waite step: subtract `round(x/2π)` copies of a
two-double 2π. That is *congruent* correct much further than it is *usable* —
the quotient is an integer, so the result stays in the right residue class
even when it is far outside `[-π, π]` — but two things degrade with `|x|`:

  * **Range.** `round(x·(1/PI2_HI))` uses a rounded reciprocal, whose ~1.1e-16
    relative error puts the quotient within 1/2 of the wrong integer once
    `|x|·1.1e-16 ≳ 1/2`. The reduced value then lands near `±(π + |x|·1.2e-16)`
    instead of inside `[-π, π]`, and Markley's starter — equation (20) is fitted
    on `M ∈ [-π, π]` — is being extrapolated.
  * **Accuracy.** 2π minus the two doubles leaves a 6.0e-33 residual, so the
    reduction carries an absolute error of about `|x|·1e-33/2π`, and above
    `|x| ≈ 2^53·2π` the quotient is no longer an exactly-representable integer
    at all and the reduction stops meaning anything.

Measured across `e ∈ {0, 0.3, 0.9}` against the scalar `Base.rem2pi` + Markley
path, the batch kernel holds the ≤3.1e-15 agreement it has at `M ≈ 1` out to
2^48, and then breaks down fast: 6.4e-15 at 2^50, 3.4e-4 at 2^52, 3.9e25 at
2^54, and `NaN` (the cubic in Markley's starter overflows) by 2^156. 2^48 is
the last magnitude with no measurable divergence at all.

This is a *domain*, not a fallback: `vrem2pi` stays branch-free, because a
range test that fell back to `Base.rem2pi` inline would put a non-vectorizable
libm call in the loop body and cost every ordinary orbit the batching. Rows
that exceed the bound are routed to the scalar solver one level up, in
`_use_simd`/`_use_dual_simd` — see `_ma_in_kernel_range`.

Note what this bound is *not*: a claim that such mean anomalies are unphysical.
`M = n·(t − tp)` exceeds 2^48 for an entirely well-posed ellipse whenever the
period is short enough or `tp` far enough from the epochs — `a ≈ 1e-11` AU at a
few thousand days, or an ordinary 1 AU orbit with `tp` 1e16 days away, both of
which a hot tempering chain reaches. `Base.rem2pi` reduces those exactly, and
after this routing so does every path in the package.
"""
const VREM2PI_MAX = 2.0^48

# rem2pi(x, RoundNearest), branch-free, fma-compensated product.
# CALLER'S CONTRACT: |x| <= VREM2PI_MAX (see above), same discipline as
# `vlog`'s positive-normal-finite domain.
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

# x^(2/3) for x ≥ 0, division-free — the Markley cubic-root step. Markley
# needs cbrt(X²); computing it as X · X^(-1/3) never materializes the square
# (no overflow ceiling at large X) and the reciprocal cube root has a
# multiply-only cubic iteration, s ← s(1 + r/3 + 2r²/9) with r = 1 − X s³
# (error → (14/3)ε³ per step; verified). Exponent-hack seed is good to 3.4%,
# so two iterations land at ≤4e-11 relative — and this only seeds E1, whose
# δ-corrections converge to the true Kepler root regardless, so the basin
# requirement is ~1e-3. Zero divides: llvm-mca shows this loop 86%
# latency-bound on register dependencies, and the previous cbrt's divides
# sat in series on the critical path.
@inline function _pow23(x::Float64)
    s = reinterpret(Float64, 0x553ef0ff289dd796 - reinterpret(UInt64, x) ÷ 3)
    r = 1 - x * s * s * s
    s = s * muladd(r, muladd(r, 2 / 9, 1 / 3), 1.0)
    r = 1 - x * s * s * s
    s = s * muladd(r, muladd(r, 2 / 9, 1 / 3), 1.0)
    return ifelse(x == 0, 0.0, x * s)
end

# Branch-free natural log for positive, normal, finite Float64 — the musl
# log.c core with the special-case branches omitted, so a loop over it
# vectorizes exactly like the kernels above. Same accuracy as the full
# routine (<1 ulp) on the domain it accepts; NaN/Inf/0/subnormal inputs are
# the *caller's* contract, exactly as `markley_sincosE` assumes e < 1.
#
# This lives here rather than in a likelihood because it is the same kind of
# object as `vsincos`: an instrument-grade math kernel whose only reason to
# exist is keeping a hot loop branch-free. First consumer is Octofitter's
# per-transit Gaussian normalization term, where the libm `log` was 13% of
# the whole DR4 primal evaluation.
#
# Only the Float64 method is the kernel; every other type (ForwardDiff Duals
# above all) falls through to `Base.log`, so a gradient evaluation never
# runs polynomial arithmetic in Dual land.
const VLN2_HI = 6.93147180369123816490e-1
const VLN2_LO = 1.90821492927058770002e-10

@inline vlog(x) = log(x)
@inline function vlog(x::Float64)
    ix = reinterpret(UInt64, x)
    # Reduce to m ∈ [√2/2, √2): subtracting the bit pattern of √2/2 makes the
    # biased exponent field count doublings from that anchor, so an arithmetic
    # shift recovers k with the mantissa's high bits acting as the rounding
    # choice between [1, 2) and [√2/2, √2) framings.
    tmp = ix - 0x3fe6a09e667f3bcd
    k = Float64(reinterpret(Int64, tmp) >> 52)
    m = reinterpret(Float64, (tmp & 0x000fffffffffffff) + 0x3fe6a09e667f3bcd)
    f = m - 1.0
    s = f / (2.0 + f)
    z = s * s
    w = z * z
    t1 = w * @evalpoly(w, 3.999999999940941908e-1, 2.222219843214978396e-1,
        1.531383769920937332e-1)
    t2 = z * @evalpoly(w, 6.666666666666735130e-1, 2.857142874366239149e-1,
        1.818357216161805012e-1, 1.479819860511658591e-1)
    hfsq = 0.5 * f * f
    return k * VLN2_HI - ((hfsq - (s * (hfsq + t2 + t1) + k * VLN2_LO)) - f)
end

# Markley (1995) without the early-return branches, on vectorizable kernels.
# Returns sincos(E) directly — the only quantities the observables need — so
# exactly ONE full sincos is evaluated per solve.
# Divides are written to minimize the number in *series*: llvm-mca on the
# vectorized loop shows 86% of backend pressure is register dependencies
# (the divide/sqrt chain), with the divider ports only ~20% loaded — so an
# algebraic rewrite that trades a divide for multiplies is a straight win
# even where it costs more instructions. Each rewrite is exact algebra or a
# same-order reformulation (identical real value up to rounding /
# truncation far below the kernel tolerance), gated by the ≤4e-15
# kernel-vs-scalar sweep in the test suite.
#
# The kernel is split at E1 into two halves so `solve_row_simd!` can run
# them as separate passes over a tile (see the fission note there); the
# fused composition below is what the Dual block path and the scalar-API
# tests use, and is the reference the halves must reproduce exactly.

# Half 1: Markley's non-iterative starter. `M` must already be reduced by
# `vrem2pi`. One divide and one sqrt.
@inline function _markley_E1(M::Float64, e::Float64)
    pi2 = π * π
    # 1/(5(1+e)) is loop-invariant (e is a row constant): as a reciprocal
    # multiply, LICM hoists it and the per-epoch divide disappears. Written
    # as x/c, the divide stays — IEEE forbids LLVM making this change.
    α = (3 * pi2 + 8 * (pi2 - π * abs(M)) * inv(5 * (1 + e))) * (1 / (pi2 - 6))
    d = 3 * (1 - e) + α * e
    q = 2 * α * d * (1 - e) - M * M
    r = 3 * α * d * (d - 1 + e) * M + M * M * M
    # abs() inside sqrt: the discriminant is provably ≥ 0 for 0 ≤ e < 1, but
    # without it LLVM keeps a domain-error throw branch that blocks the
    # loop vectorizer.
    w = _pow23(abs(r) + sqrt(abs(q * q * q + r * r)))
    # (2rw/den + M)/d = (2rw + M·den)/(d·den): two serial divides → one.
    den = @evalpoly(w, q * q, q, 1)
    return (2 * r * w + M * den) / (d * den)
end

# Half 2: the 4th-order correction, as a series reversion of Kepler's
# equation about E1 rather than the nested δ3→δ4→δ5 evaluation. The nested
# form costs three divides in series; the reverted series costs ONE
# reciprocal (1/f1) with everything else multiply-adds:
#
#     δ = −u(1 + Au + (2A² − B)u² + (5A³ − 5AB + C)u³),
#     u = f0/f1,  A = f2/2f1,  B = f3/6f1,  C = f4/24f1 = −A/12
#
# (f4 = −f2 for Kepler's equation). Same formal order as δ5; the truncated
# cross terms are O(u⁵) with u ≲ 1e-3 from Markley's starter, far below the
# kernel tolerance.
@inline function _markley_correct(M::Float64, E1::Float64, e::Float64)
    sE1, cE1 = vsincos(E1)
    f2 = e * sE1
    f3 = e * cE1
    f0 = E1 - f2 - M
    f1 = 1 - f3
    invf1 = inv(f1)
    u = f0 * invf1
    A = 0.5 * f2 * invf1
    B = f3 * invf1 * (1 / 6)
    A2 = A * A
    δ = -u * @evalpoly(u, 1.0, A, 2 * A2 - B, 5 * A * A2 - 5 * A * B - A * (1 / 12))
    # sincos(E1 + δ) via small-angle rotation — no second full sincos
    sδ, cδ = vsincos_small(δ)
    sE = sE1 * cδ + cE1 * sδ
    cE = cE1 * cδ - sE1 * sδ
    return sE, cE
end

@inline function markley_sincosE(Min::Float64, e::Float64)
    M = vrem2pi(Min)
    E1 = _markley_E1(M, e)
    return _markley_correct(M, E1, e)
end

# Epochs per fission tile. The fused loop body was ~240 µops — the ROB
# holds barely two iterations, which is why llvm-mca measured 86%
# register-dependency stalls with the divider ports mostly idle. Splitting
# the kernel at E1 halves each body, so several iterations of each pass
# overlap in flight; the seed values live in two tile-sized stack (MVector)
# buffers, costing no trajectory storage. MUST stay a multiple of
# `DUAL_SIMD_BLOCK`: chunk boundaries in the threaded solve are aligned to
# that block, and a tile that respects it leaves the vector-body/scalar-tail
# split of every epoch invariant across serial and chunked runs — which is
# the "bit for bit identical to serial" contract in `orbitsolve!`.
const SIMD_TILE = 32

# SIMD Float64 batch path for one hierarchy row: pass 1 seeds (rem2pi +
# Markley starter), pass 2 corrects and stores states.
function solve_row_simd!(traj::Trajectory{Float64}, row::Row{Float64}, j::Int)
    n_per_day = row.n / year2day_julian
    tp = row.tp
    t_em = traj.t_em
    rx = traj.rx; ry = traj.ry; rz = traj.rz
    rvx = traj.rvx; rvy = traj.rvy; rvz = traj.rvz
    e = row.e
    nep = length(t_em)
    Mb = MVector{SIMD_TILE,Float64}(undef)
    Eb = MVector{SIMD_TILE,Float64}(undef)
    k0 = 0
    @inbounds while k0 < nep
        len = min(SIMD_TILE, nep - k0)
        @simd ivdep for i in 1:len
            M = vrem2pi(n_per_day * (t_em[k0+i] - tp))
            Mb[i] = M
            Eb[i] = _markley_E1(M, e)
        end
        @simd ivdep for i in 1:len
            sE, cE = _markley_correct(Mb[i], Eb[i], e)
            x, y, z, vx, vy, vz = _states_from_E(row, sE, cE)
            k = k0 + i
            rx[k, j] = x; ry[k, j] = y; rz[k, j] = z
            rvx[k, j] = vx; rvy[k, j] = vy; rvz[k, j] = vz
        end
        k0 += len
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
