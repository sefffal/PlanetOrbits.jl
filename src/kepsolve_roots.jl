using .Roots: Roots

@inline function kepler_solver(_MA::Real, e::Real, method::RootsMethod)
    MA = rem2pi(_MA, RoundNearest)

    # Standard Kepler's equation and derivatives
    kep1(EA) = EA - MA - e*sin(EA)
    kep1′(EA) = 1 - e*cos(EA)
    kep1′′(EA) = e*sin(EA)
    kep1′′′(EA) = e*cos(EA)
    kep1′′′′(EA) = -e*sin(EA)
    kep1′′′′′(EA) = -e*cos(EA)

    fs = (kep1,kep1′,kep1′′,kep1′′′,kep1′′′′,kep1′′′′)

    # # Hyperbolic Kepler's equation
    # kep2(EA) = e*sinh(EA) - EA - MA
    # if e < 1
    #     f = kep1
    # else
    #     f = kep2
    # end
    if e == 0
        return MA
    end

    if typeof(method.method) <: Roots.AbstractBracketingMethod
        # if -π < MA < 0 || π < MA
        #     initial = 
        # else
        #     initial = MA + e
        # end
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

