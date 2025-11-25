
"""
    PlanetOrbits.Markley()

Kepler solver implementation from AstroLib, based on Markley (1995)
Celestial Mechanics and Dynamical Astronomy, 63, 101 (DOI:10.1007/BF00691917). 

"""
struct Markley <: AbstractSolver end

# Our own simpler version of rem2pi. This probably doesn't handle
# inf/nan and the boundary conditions perfectly.
# However, it works with Enzyme
function myrem2pi(x)                                                                                                                                                
    if x > π                                                                                                                                                        
        while x > π                                                                                                                                                 
            x -= 2π                                                                                                                                                 
        end                                                                                                                                                         
    elseif x < -π                                                                                                                                                   
        while x < -π                                                                                                                                                
            x += 2π                                                                                                                                                 
        end                                                                                                                                                         
    end                                                                                                                                                             
    x                                                                                                                                                               
end

# The following function is taken directly from AstroLib.jl
# We remove one invariant check we handle elsewhere and also
# force inlining for about a 5% speedup.
# We also supply analytic gradients for use in autodiff packages.
@inline function kepler_solver(_M::Real, e::Real, ::Markley)
    # We already handle this invariant
    # @assert 0 <= e <= 1 "eccentricity must be in the range [0, 1]"
    # M must be in the range [-pi, pi], see Markley (1995), page 2.
    M = rem2pi(_M, RoundNearest)
    # M = myrem2pi(_M)
    T = float(promote_type(typeof(M), typeof(e)))
    if iszero(M) || iszero(e)
        return T(M)
    end
    pi2 = abs2(T(pi))
    # equation (20)
    α = (3 * pi2 + 8 * (pi2 - pi * abs(M)) / (5 * (1 + e)))/(pi2 - 6)
    # equation (5)
    d = 3 * (1 - e) + α * e
    # equation (9)
    q = 2 * α * d * (1 - e) - M * M
    # equation (10)
    r = 3 * α * d * (d - 1 + e) * M + M * M * M
    # equation (14)
    w = cbrt(abs2(abs(r) + sqrt(q * q * q + r * r)))
    # equation (15)
    E1 = (2 * r * w / @evalpoly(w, q * q, q, 1) + M)/d
    # equation (26) & equation (27)
    f2, f3 = e .* sincos(E1)
    # equation (21)
    f0 = E1 - f2 - M
    # equation (25)
    f1 = 1 - f3
    # equation (22)
    δ3 = -f0 / (f1 - f0 * f2 / (2 * f1))
    # equation (23)
    δ4 = -f0 / @evalpoly(δ3, f1, f2 / 2, f3 / 6)
    # equations (24) and (28)
    δ5 = -f0 / @evalpoly(δ4, f1, f2 / 2, f3 / 6, - f2 / 24)
    return E1 + δ5 # equation 29
end

"""
    kepler_solver_phased!(EA, MA, e, M_buf, E1_buf)

Vectorized Kepler solver using the Markley (1995) algorithm for multiple mean anomalies.
Operates in-place on the output vector EA.

This implementation separates vectorizable and non-vectorizable operations into phases
to maximize SIMD utilization:
1. Phase 1: Range reduction (fully vectorizable)
2. Phase 2: Algebraic computation through cbrt (partially vectorizable)
3. Phase 3: Trigonometric refinement (vectorizable with SIMD)

# Arguments
- `EA::AbstractVector{T}`: Output vector for eccentric anomalies (modified in-place)
- `MA::AbstractVector{T}`: Input vector of mean anomalies
- `e::T`: Eccentricity (scalar, same for all anomalies)
- `M_buf::AbstractVector{T}`: Workspace buffer for range-reduced mean anomalies
- `E1_buf::AbstractVector{T}`: Workspace buffer for initial eccentric anomaly estimates

All vectors must have the same length. Workspace buffers can be reused across calls.

# Returns
The modified `EA` vector containing the eccentric anomalies.

# Example
```julia
MA = [0.0, 1.0, 2.0, 3.0]
e = 0.5
EA = similar(MA)
M_buf = similar(MA)
E1_buf = similar(MA)
kepler_solver_phased!(EA, MA, e, M_buf, E1_buf)
```
"""
function kepler_solver_phased!(EA::AbstractVector{T}, MA::AbstractVector{T}, e::T,
                            M_buf::AbstractVector{T}, E1_buf::AbstractVector{T}) where {T<:AbstractFloat}
    n = length(MA)
    @assert length(EA) == n == length(M_buf) == length(E1_buf) "All vectors must have the same length"

    inv2π = T(inv(2π))
    twopi = T(2π)

    # Fast path for circular orbits
    if abs(e) < eps(T)
        @inbounds @simd for i in 1:n
            M_raw = MA[i]
            EA[i] = M_raw - twopi * round(M_raw * inv2π)
        end
        return EA
    end

    # Precompute eccentricity-dependent constants
    π_T = T(π)
    π2 = π_T * π_T
    denom_α = π2 - T(6)
    one_plus_e = one(T) + e
    one_minus_e = one(T) - e
    d0 = T(3) * one_minus_e
    coeff_α1 = T(3) * π2
    coeff_α2 = T(8) / (T(5) * one_plus_e)

    two = T(2)
    three = T(3)
    half = T(0.5)
    sixth = T(1)/T(6)
    inv24 = T(1)/T(24)

    # Phase 1: Range reduction - fully SIMD
    @inbounds @simd for i in 1:n
        M_raw = MA[i]
        M_buf[i] = M_raw - twopi * round(M_raw * inv2π)
    end

    # Phase 2: Compute E1 through algebraic formula
    @inbounds @fastmath @simd for i in 1:n
        M = M_buf[i]
        absM = abs(M)
        α = (coeff_α1 + coeff_α2 * (π2 - π_T * absM)) / denom_α
        d = d0 + α * e
        q = two * α * d * one_minus_e - M * M
        r = three * α * d * (d - one_minus_e) * M + M * M * M

        abs_r = abs(r)
        discriminant = q * q * q + r * r
        sqrt_disc = sqrt(abs(discriminant))
        w_inner = abs_r + sqrt_disc
        w = cbrt(w_inner * w_inner)  # x^(2/3) directly
        denom_E1 = w*w + q*w + q*q
        E1_buf[i] = (two * r * w / denom_E1 + M) / d
    end

    # Phase 3: Newton refinement with trig functions
    @inbounds @fastmath @simd for i in 1:n
        M = M_buf[i]
        E1 = E1_buf[i]

        sinE1 = sin(E1)
        cosE1 = cos(E1)
        f2 = e * sinE1
        f3 = e * cosE1
        f0 = E1 - f2 - M
        f1 = one(T) - f3

        δ3 = -f0 / (f1 - f0 * f2 / (two * f1))
        δ4 = -f0 / (f1 + δ3 * f2 * half + δ3 * δ3 * f3 * sixth)
        δ4_sq = δ4 * δ4
        δ5 = -f0 / (f1 + δ4 * f2 * half + δ4_sq * f3 * sixth - δ4_sq * δ4 * f2 * inv24)

        EA[i] = E1 + δ5
    end

    return EA
end

