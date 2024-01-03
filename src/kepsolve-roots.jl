using Roots

@inline function PlanetOrbits.kepler_solver(_MA::Real, e::Real, method::RootsMethod)
    if e < 1
        MA = rem2pi(_MA, RoundNearest)
        return kepler_solver_roots(MA, e, method)
    else
        @show _MA
        return hyperbolic_kepler_solver_roots(_MA, e, method)
    end
end

@inline function kepler_solver_roots(MA::Real, e::Real, method::RootsMethod)

    # Standard Kepler's equation and derivatives
    kep1(EA) = EA - MA - e*sin(EA)
    kep1′(EA) = 1 - e*cos(EA)
    kep1′′(EA) = e*sin(EA)
    kep1′′′(EA) = e*cos(EA)
    kep1′′′′(EA) = -e*sin(EA)
    kep1′′′′′(EA) = -e*cos(EA)

    fs = (kep1,kep1′,kep1′′,kep1′′′,kep1′′′′,kep1′′′′)

    if e == 0
        return MA
    end

    if typeof(method.method) <: Roots.AbstractBracketingMethod
        if e < 0.7
            initial = (MA - e, MA+e)
        else
            initial = machin(e, MA) .+ (-0.05, 0.05)
        end
    else
        if e < 0.5
            initial = MA
        else
            if -π < MA < 0 || π < MA
                initial = MA - e
            else
                initial = MA + e
            end
        end
    end

    EA = Roots.find_zero(fs, initial, method.method; method.kwargs...)
    
    return EA
end


@inline function hyperbolic_kepler_solver_roots(MA::Real, e::Real, method::RootsMethod)

    # Hyperbolic Kepler's equation and derivatives
    keph(EA) = - EA - MA + e*sinh(EA)
    keph′(EA) = e*cosh(EA) - 1
    keph′′(EA) = e*sinh(EA)
    keph′′′(EA) = e*cosh(EA)
    keph′′′′(EA) = e*sinh(EA)
    keph′′′′′(EA) = e*cosh(EA)

    fs = (keph,keph′,keph′′,keph′′′,keph′′′′,keph′′′′)


    if typeof(method.method) <: Roots.AbstractBracketingMethod
        initial = (MA - e, MA+e)
    else
        if -π < MA < 0 || π < MA
            initial = MA - e
        else
            initial = MA + e
        end
    end

    EA = Roots.find_zero(fs, initial, method.method; method.kwargs...)
    
    return EA
end

# Functions to find a good starting point for iteration,
# specifically with Newton's method.
# Functions by J Cook: https://www.johndcook.com/blog/2022/11/02/keplers-equation-python/

# These are currently unused. For Newton's method, we are faster
# just running more iterations from a simple guess than calculating
# all these sqrt and cbrt to get a good starting point.

# This will solve the special form of the cubic we need.
function solve_cubic(a, c, d)
    p = c/a
    q = d/a
    k = sqrt( q^2/4 + p^3/27 )
    return cbrt(-q/2 - k) + cbrt(-q/2 + k)
end

# Machin's starting point for Newton's method
# See johndcook.com/blog/2022/11/01/kepler-newton/
function machin(e, M)
    n = sqrt(5 + sqrt(16 + 9/e))
    a = n*(e*(n^2 - 1)+1)/6
    c = n*(1-e)
    d = -M
    s = solve_cubic(a, c, d)
    return n*asin(s)  
end  
