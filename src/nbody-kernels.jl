# ---------------------------------------------------
# N-body kernels: the AHL21 symplectic map on StaticArrays state.
#
# Derived from NbodyGradient.jl v0.2.3 (7d9225b), © 2021 Eric Agol, David
# Hernandez & Zach Langford, MIT license — merged into PlanetOrbits with the
# authors' agreement. The numerical method (universal-variable Kepler drift,
# G/H series, compensated summation, map composition) is preserved verbatim;
# only the state representation differs: instead of a fat mutable `State`,
# the integrator state is an immutable `NBodyState` of `SMatrix` values and
# every step is a pure `state -> state` function. Kernel equivalence against
# upstream is enforced by checked-in fixtures (test/fixtures/nbody_reference.jl).
#
# Algorithm: Agol, Hernandez & Langford (2021), MNRAS 507, 1582
# (arXiv:2106.02188). Please cite this paper when publishing results computed
# with the AHL21 propagator.
#
# Conventions here: G is *parametrized* — every kernel takes the per-body
# gravitational-mass vector `Gm` (= G .* masses in whatever unit system the
# caller uses; the AHL21 propagator uses AU, days, M⊙). Upstream's kick-group
# machinery (`pair`) is not carried over: every pair is treated as a
# Keplerian pair, which is upstream's default configuration (all-false
# `pair`) and makes upstream's `kickfast!`/`phic!` no-ops.
# ---------------------------------------------------

"""
    NBodyState{N,T}

Immutable integrator state for `N` bodies: barycentric positions [AU],
velocities [AU/day], and the Kahan compensated-summation residuals carried
between steps. A pure value — stepping produces a new state.
"""
struct NBodyState{N,T,L}
    x::SMatrix{3,N,T,L}
    v::SMatrix{3,N,T,L}
    xerr::SMatrix{3,N,T,L}
    verr::SMatrix{3,N,T,L}
end

NBodyState(x::SMatrix{3,N,T,L}, v::SMatrix{3,N,T,L}) where {N,T,L} =
    NBodyState{N,T,L}(x, v, zero(x), zero(v))

# Time reversal: marching backwards is implemented as forward marching of the
# velocity-reversed state (the map is time-reversal equivariant), keeping the
# universal-variable γ positive.
_reverse_v(st::NBodyState{N,T,L}) where {N,T,L} =
    NBodyState{N,T,L}(st.x, -st.v, st.xerr, -st.verr)

"""
    ahl21_step(st::NBodyState, Gm::SVector, h) -> NBodyState

One step of the AHL21 symplectic map (Agol, Hernandez & Langford 2021):

    drift(h/2) ∘ kepler-drift pairs(h/2, forward order) ∘ ϕ†(h)
              ∘ kepler-drift pairs(h/2, reverse order) ∘ drift(h/2)

with compensated summation throughout. Second-order, symplectic, exact for a
two-body system (the fourth-order corrector ϕ† vanishes identically for one
pair). `Gm` is G×mass per body, in units consistent with `st` and `h`.
"""
function ahl21_step(st::NBodyState{N,T,L}, Gm::SVector{N}, h) where {N,T,L}
    x = MMatrix(st.x); xerr = MMatrix(st.xerr)
    v = MMatrix(st.v); verr = MMatrix(st.verr)
    h2 = 0.5 * h
    _drift!(x, xerr, v, h2)
    _drift_kepler!(x, v, xerr, verr, Gm, h2)
    _phisalpha!(x, v, verr, Gm, h)
    _kepler_drift!(x, v, xerr, verr, Gm, h2)
    _drift!(x, xerr, v, h2)
    return NBodyState{N,T,L}(SMatrix(x), SMatrix(v), SMatrix(xerr), SMatrix(verr))
end

# Kahan (1965) compensated summation: returns the updated (sum, residual).
@inline function _comp_sum(sum_value, sum_error, addend)
    sum_error += addend
    tmp = sum_value + sum_error
    sum_error += sum_value - tmp
    return tmp, sum_error
end

@inline _dot3(a, b) = a[1]*b[1] + a[2]*b[2] + a[3]*b[3]
@inline _dot3(a) = _dot3(a, a)

# Free drift of every body.
@inline function _drift!(x, xerr, v, h)
    @inbounds for i in axes(x, 2), k in 1:3
        x[k, i], xerr[k, i] = _comp_sum(x[k, i], xerr[k, i], h * v[k, i])
    end
    return nothing
end

# Combined -drift+Kepler over all pairs, forward order (start of the step).
@inline function _drift_kepler!(x, v, xerr, verr, Gm::SVector{N}, h) where {N}
    @inbounds for i in 1:N-1, j in i+1:N
        _kepler_driftij!(x, v, xerr, verr, Gm, i, j, h, true)
    end
    return nothing
end

# Combined Kepler+(-drift) over all pairs, reverse order (end of the step).
@inline function _kepler_drift!(x, v, xerr, verr, Gm::SVector{N}, h) where {N}
    @inbounds for i in N-1:-1:1, j in N:-1:i+1
        _kepler_driftij!(x, v, xerr, verr, Gm, i, j, h, false)
    end
    return nothing
end

# Kepler step minus linear drift for the pair (i, j), distributed onto the
# two bodies by mass fraction (upstream `kepler_driftij_gamma!`).
@inline function _kepler_driftij!(x, v, xerr, verr, Gm, i::Int, j::Int, h, drift_first::Bool)
    @inbounds x0 = SVector(x[1,i]-x[1,j], x[2,i]-x[2,j], x[3,i]-x[3,j])
    @inbounds v0 = SVector(v[1,i]-v[1,j], v[2,i]-v[2,j], v[3,i]-v[3,j])
    gm = Gm[i] + Gm[j]
    iszero(gm) && return nothing   # two test particles: pure drift, no correction
    delx, delv = _delxv_gamma(x0, v0, gm, h, drift_first)
    gminv = inv(gm)
    mi = Gm[i] * gminv             # mass fractions (G cancels)
    mj = Gm[j] * gminv
    @inbounds for k in 1:3
        x[k,i], xerr[k,i] = _comp_sum(x[k,i], xerr[k,i],  mj*delx[k])
        x[k,j], xerr[k,j] = _comp_sum(x[k,j], xerr[k,j], -mi*delx[k])
        v[k,i], verr[k,i] = _comp_sum(v[k,i], verr[k,i],  mj*delv[k])
        v[k,j], verr[k,j] = _comp_sum(v[k,j], verr[k,j], -mi*delv[k])
    end
    return nothing
end

# The fourth-order corrector ϕ† ("phisalpha" with α = 2), applied to every
# pair. Vanishes identically for a single pair (two-body systems).
# (@inline like every MMatrix-mutating kernel: they must all inline into
# ahl21_step so its state buffers provably don't escape and stay on-stack.)
@inline function _phisalpha!(x::MMatrix{3,N,T}, v, verr, Gm, h) where {N,T}
    a = zero(MMatrix{3,N,T})
    @inbounds for i in 1:N-1, j in i+1:N
        rij = SVector(x[1,i]-x[1,j], x[2,i]-x[2,j], x[3,i]-x[3,j])
        r2 = _dot3(rij)
        r3 = r2 * sqrt(r2)
        for k in 1:3
            fac = rij[k] / r3
            a[k,i] -= Gm[j] * fac
            a[k,j] += Gm[i] * fac
        end
    end
    coeff = h^3 / 24    # upstream α h³/96 · 2 with α = 2, G folded into Gm
    @inbounds for i in 1:N-1, j in i+1:N
        aij = SVector(a[1,i]-a[1,j], a[2,i]-a[2,j], a[3,i]-a[3,j])
        rij = SVector(x[1,i]-x[1,j], x[2,i]-x[2,j], x[3,i]-x[3,j])
        r2 = _dot3(rij)
        r1 = sqrt(r2)
        ardot = _dot3(aij, rij)
        fac1 = coeff / (r2 * r2 * r1)
        fac2 = 2 * (Gm[i] + Gm[j]) / r1 + 3 * ardot
        for k in 1:3
            fac = fac1 * (rij[k] * fac2 - r2 * aij[k])
            v[k,i], verr[k,i] = _comp_sum(v[k,i], verr[k,i],  Gm[j]*fac)
            v[k,j], verr[k,j] = _comp_sum(v[k,j], verr[k,j], -Gm[i]*fac)
        end
    end
    return nothing
end

# ---------------------------------------------------
# Universal-variable Kepler drift (upstream `jac_delxv_gamma!`, no-gradient
# path): evolve the relative state (x0, v0) of one pair under two-body
# gravity with GM = gm for time h, minus the linear drift, expressed through
# Gauss f/g functions of the universal anomaly γ. `drift_first` selects
# whether the removed drift precedes (true) or follows (false) the Kepler
# step. Returns (Δx, Δv).
# ---------------------------------------------------
function _delxv_gamma(x0::SVector{3}, v0::SVector{3}, gm, h, drift_first::Bool)
    rtmp = x0 - h * v0
    r0 = drift_first ? sqrt(_dot3(rtmp)) : sqrt(_dot3(x0))
    r0inv = inv(r0)
    beta0 = 2 * gm * r0inv - _dot3(v0)
    beta0inv = inv(beta0)
    signb = sign(beta0)
    sqb = sqrt(signb * beta0)
    zeta = gm - r0 * beta0
    eta = drift_first ? _dot3(rtmp, v0) : _dot3(x0, v0)
    # Initial guess for γ (cubic/quadratic in the small-γ expansion):
    if !iszero(zeta)
        zinv = 6 / zeta
        gamma = _cubic1(0.5 * eta * sqb * zinv, r0 * signb * beta0 * zinv,
                        -h * signb * beta0 * sqb * zinv)
    elseif !iszero(eta)
        reta = r0 / eta
        disc = reta^2 + 2h / eta
        gamma = disc > 0 ? sqb * (-reta + sqrt(disc)) : h * r0inv * sqb
    else
        gamma = h * r0inv * sqb
    end
    # Newton solve of the universal Kepler equation for γ:
    gamma1 = 2gamma
    gamma2 = 3gamma
    c2 = -2zeta
    c3 = 2 * eta * signb * sqb
    c4 = -sqb * h * beta0
    iter = 0
    while iter < 20
        gamma2 = gamma1
        gamma1 = gamma
        xx = 0.5 * gamma
        if beta0 > 0
            sx, cx = sincos(xx)
        else
            sx = sinh(xx); cx = exp(-xx) + sx
        end
        gamma -= (gm*gamma + c2*sx*cx + c3*sx^2 + c4) /
                 (2*signb*zeta*sx^2 + c3*sx*cx + r0*beta0)
        iter += 1
        (gamma == gamma2 || gamma == gamma1) && break
    end
    # Wisdom/Hernandez Gᵢ(β, γ) functions at the converged γ:
    xx = 0.5 * gamma
    if beta0 > 0
        sx, cx = sincos(xx)
    else
        sx = sinh(xx); cx = exp(-xx) + sx
    end
    g1bs = 2 * sx * cx / sqb
    g2bs = 2 * signb * sx^2 * beta0inv
    g0bs = 1 - beta0 * g2bs
    g3bs = _G3(gamma, beta0, sqb)
    r = r0 * g0bs + eta * g1bs + gm * g2bs
    rinv = inv(r)
    dfdt = -gm * g1bs * rinv * r0inv
    if drift_first
        fm1 = -gm * r0inv * g2bs                        # f - 1
        gmh = gm * r0inv * (h * g2bs - r0 * g3bs)       # g - h f
        dgdtm1 = gm * r0inv * rinv * (h * g1bs - r0 * g2bs)  # ġ - h ḟ - 1
    else
        h1 = _H1(gamma, beta0)
        h2 = _H2(gamma, beta0, sqb)
        fm1 = gm * rinv * (g2bs - gm * r0inv * h1)      # f - 1 - h ḟ
        gmh = gm * rinv * (r0 * h2 + eta * h1)          # g - h ġ
        dgdtm1 = -gm * rinv * g2bs                      # ġ - 1
    end
    delx = fm1 * x0 + gmh * v0
    delv = dfdt * x0 + dgdtm1 * v0
    return delx, delv
end

# Real root of x³ + a x² + b x + c (upstream `cubic1`), used for the γ guess.
function _cubic1(a, b, c)
    third = one(a) / 3
    a3 = a * third
    Q = a3^2 - b * third
    R = a3^3 + 0.5 * (-a3 * b + c)
    R2 = R^2
    Q3 = Q^3
    if R2 < Q3
        return -c / b
    else
        A = -sign(R) * cbrt(abs(R) + sqrt(R2 - Q3))
        B = iszero(A) ? zero(A) : Q / A
        return A + B - a3
    end
end

# ---------------------------------------------------
# G₃/H₁/H₂ special functions with small-γ series (upstream `G3`,`H1`,`H2`).
# The series are iterated to fixed point in the working precision, matching
# upstream exactly (including the termination conditions).
# ---------------------------------------------------
function _G3(gamma, beta, sqb; gc=oftype(gamma, 0.5))
    if gamma < gc
        return _G3_series(gamma, beta, sqb)
    elseif beta >= 0
        return (gamma - sin(gamma)) / (sqb * beta)
    else
        return (gamma - sinh(gamma)) / (sqb * beta)
    end
end

function _G3_series(gamma, beta, sqb)
    x2 = -sign(beta) * gamma^2
    term = one(gamma)
    g3 = one(gamma)
    g31 = 2g3
    g32 = 2g3
    n = 0
    iter = 0
    while true
        g32 = g31
        g31 = g3
        n += 1
        term *= x2 / ((2n + 3) * (2n + 2))
        g3 += term
        iter += 1
        (iter >= 100 || g3 == g32 || g3 == g31) && break
    end
    return g3 * (-x2 * gamma / (6 * beta * sqb))
end

function _H1(gamma, beta; gc=oftype(gamma, 0.5))
    if gamma < gc
        return _H1_series(gamma, beta)
    elseif beta >= 0
        return (4 * sin(0.5 * gamma)^2 - gamma * sin(gamma)) / beta^2
    else
        return (-4 * sinh(0.5 * gamma)^2 + gamma * sinh(gamma)) / beta^2
    end
end

function _H1_series(gamma, beta)
    x2 = -sign(beta) * gamma^2
    term = one(gamma)
    h1 = one(gamma)
    h11 = 2h1
    h12 = 2h1
    n = 0
    iter = 0
    while true
        h12 = h11
        h11 = h1
        n += 1
        term *= x2 * (n + 1)
        term /= (2n + 4) * (2n + 3) * n
        h1 += term
        iter += 1
        (iter >= 100 || h1 == h12 || h1 == h11) && break
    end
    return h1 * (x2^2 / (12 * beta^2))
end

function _H2(gamma, beta, sqb; gc=oftype(gamma, 0.5))
    if gamma < gc
        return _H2_series(gamma, beta, sqb)
    elseif beta >= 0
        return (sin(gamma) - gamma * cos(gamma)) / (sqb * beta)
    else
        return (sinh(gamma) - gamma * cosh(gamma)) / (sqb * beta)
    end
end

function _H2_series(gamma, beta, sqb)
    x2 = -sign(beta) * gamma^2
    term = one(gamma)
    h2 = one(gamma)
    h21 = 2h2
    h22 = 2h2
    n = 0
    iter = 0
    while true
        h22 = h21
        h21 = h2
        n += 1
        term *= x2
        term /= (4n + 6) * n
        h2 += term
        iter += 1
        (iter >= 100 || h2 == h22 || h2 == h21) && break
    end
    return h2 * (-x2 * gamma / (3 * beta * sqb))
end
