
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
    # M = rem2pi(_M, RoundNearest)
    M = myrem2pi(_M)
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

