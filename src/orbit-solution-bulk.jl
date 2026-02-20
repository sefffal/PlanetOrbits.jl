"""
    OrbitSolutionBulk{T, TArr, TEl, TTime}

Bulk orbit solution storing per-epoch arrays of computed orbital quantities.
`elem` is a single orbit object (scalar fields), while ν, sinν_ω, etc. are
arrays of length N_epochs.

This is the array-of-structs → struct-of-arrays transpose of
`Vector{OrbitSolutionKep}`, designed for Reactant/XLA traceability.
"""
struct OrbitSolutionBulk{T, TArr<:AbstractVector{T}, TEl<:AbstractOrbit, TTime<:AbstractVector}
    elem::TEl
    ν::TArr
    sinν_ω::TArr
    cosν_ω::TArr
    ecosν::TArr
    r::TArr
    t::TTime
end

"""
    orbitsolve_bulk(elem::KepOrbit, times::AbstractVector)

Solve a Keplerian orbit at multiple epochs, returning an `OrbitSolutionBulk`.
Only supports elliptic orbits (e < 1). Uses the branchless Markley Kepler
solver for Reactant/XLA traceability.
"""
function orbitsolve_bulk(elem::KepOrbit{Te}, times::AbstractVector) where Te
    e = eccentricity(elem)
    T = float(promote_type(Te, eltype(times)))
    N = length(times)

    ν      = Vector{T}(undef, N)
    sinν_ω = Vector{T}(undef, N)
    cosν_ω = Vector{T}(undef, N)
    ecosν  = Vector{T}(undef, N)
    r      = Vector{T}(undef, N)

    tₚ = periastron(elem)
    n  = meanmotion(elem)

    for j in eachindex(times)
        # Mean anomaly
        MA = n / year2day_julian * (times[j] - tₚ)

        # Eccentric anomaly via branchless Markley
        EA = kepler_solver_reactant(MA, e)

        # True anomaly via half-angle formula (elliptic only, branchless)
        β = e / (1 + sqrt(1 - e^2))
        sea, cea = sincos(EA)
        ν[j] = EA + 2 * atan(β * sea, 1 - β * cea)

        # Quantities needed by downstream observables
        sinν_ω[j], cosν_ω[j] = sincos(elem.ω + ν[j])
        ecosν[j] = e * cos(ν[j])
        r[j] = elem.p / (1 + ecosν[j])
    end

    return OrbitSolutionBulk(elem, ν, sinν_ω, cosν_ω, ecosν, r, times)
end

"""
    orbitsolve_bulk(elem::RadialVelocityOrbit, times::AbstractVector)

Solve a radial-velocity orbit at multiple epochs, returning an `OrbitSolutionBulk`.
Only supports elliptic orbits (e < 1). Uses the branchless Markley Kepler
solver for Reactant/XLA traceability.
"""
function orbitsolve_bulk(elem::RadialVelocityOrbit{Te}, times::AbstractVector) where Te
    e = eccentricity(elem)
    T = float(promote_type(Te, eltype(times)))
    N = length(times)

    ν      = Vector{T}(undef, N)
    sinν_ω = Vector{T}(undef, N)
    cosν_ω = Vector{T}(undef, N)
    ecosν  = Vector{T}(undef, N)
    r      = Vector{T}(undef, N)

    tₚ = periastron(elem)
    n  = meanmotion(elem)
    p  = elem.a * (1 - e^2)

    for j in eachindex(times)
        MA = n / year2day_julian * (times[j] - tₚ)
        EA = kepler_solver_reactant(MA, e)

        β = e / (1 + sqrt(1 - e^2))
        sea, cea = sincos(EA)
        ν[j] = EA + 2 * atan(β * sea, 1 - β * cea)

        sinν_ω[j], cosν_ω[j] = sincos(elem.ω + ν[j])
        ecosν[j] = e * cos(ν[j])
        r[j] = p / (1 + ecosν[j])
    end

    return OrbitSolutionBulk(elem, ν, sinν_ω, cosν_ω, ecosν, r, times)
end

# Delegate to parent orbit for Visual wrappers
function orbitsolve_bulk(vis::VisualOrbit, times::AbstractVector)
    return orbitsolve_bulk(vis.parent, times)
end


# --- Observable functions for OrbitSolutionBulk ---

# Radial velocity [m/s]
function radvel(sol::OrbitSolutionBulk)
    K = sol.elem.K
    ecosω = sol.elem.ecosω
    rv = similar(sol.cosν_ω)
    for j in eachindex(rv)
        rv[j] = K * (sol.cosν_ω[j] + ecosω)
    end
    return rv
end

# Position in x direction [AU]
function posx(sol::OrbitSolutionBulk)
    sinΩ = sol.elem.sinΩ
    cosi_cosΩ = sol.elem.cosi_cosΩ
    x = similar(sol.r)
    for j in eachindex(x)
        x[j] = sol.r[j] * (sol.cosν_ω[j] * sinΩ + sol.sinν_ω[j] * cosi_cosΩ)
    end
    return x
end

# Position in y direction [AU]
function posy(sol::OrbitSolutionBulk)
    cosΩ = sol.elem.cosΩ
    cosi_sinΩ = sol.elem.cosi_sinΩ
    y = similar(sol.r)
    for j in eachindex(y)
        y[j] = sol.r[j] * (sol.cosν_ω[j] * cosΩ - sol.sinν_ω[j] * cosi_sinΩ)
    end
    return y
end

# Position in z direction [AU]
function posz(sol::OrbitSolutionBulk)
    sini = sol.elem.sini
    z = similar(sol.r)
    for j in eachindex(z)
        z[j] = sol.r[j] * sol.sinν_ω[j] * sini
    end
    return z
end

# Velocity in x direction [AU/julian year]
function velx(sol::OrbitSolutionBulk)
    J = sol.elem.J
    cosi_cosΩ = sol.elem.cosi_cosΩ
    sinΩ = sol.elem.sinΩ
    ecosω = sol.elem.ecosω
    esinω = sol.elem.esinω
    result = similar(sol.cosν_ω)
    for j in eachindex(result)
        result[j] = J * (cosi_cosΩ * (sol.cosν_ω[j] + ecosω) - sinΩ * (sol.sinν_ω[j] + esinω))
    end
    return result
end

# Velocity in y direction [AU/julian year]
function vely(sol::OrbitSolutionBulk)
    J = sol.elem.J
    cosi_sinΩ = sol.elem.cosi_sinΩ
    cosΩ = sol.elem.cosΩ
    ecosω = sol.elem.ecosω
    esinω = sol.elem.esinω
    result = similar(sol.cosν_ω)
    for j in eachindex(result)
        result[j] = -J * (cosi_sinΩ * (sol.cosν_ω[j] + ecosω) + cosΩ * (sol.sinν_ω[j] + esinω))
    end
    return result
end

# Velocity in z direction [AU/julian year]
function velz(sol::OrbitSolutionBulk)
    rv = radvel(sol)
    for j in eachindex(rv)
        rv[j] *= m2au * year2sec_julian
    end
    return rv
end


# --- Mass-ratio methods (star's motion due to planet) ---
# These compute the primary's observable = -M_planet/M_tot * secondary's observable

function radvel(sol::OrbitSolutionBulk, M_planet)
    rv = radvel(sol)
    M_tot = totalmass(sol.elem)
    scale = -M_planet / M_tot
    for j in eachindex(rv)
        rv[j] *= scale
    end
    return rv
end

function posx(sol::OrbitSolutionBulk, M_planet)
    x = posx(sol)
    M_tot = totalmass(sol.elem)
    scale = -M_planet / M_tot
    for j in eachindex(x)
        x[j] *= scale
    end
    return x
end

function posy(sol::OrbitSolutionBulk, M_planet)
    y = posy(sol)
    M_tot = totalmass(sol.elem)
    scale = -M_planet / M_tot
    for j in eachindex(y)
        y[j] *= scale
    end
    return y
end

function posz(sol::OrbitSolutionBulk, M_planet)
    z = posz(sol)
    M_tot = totalmass(sol.elem)
    scale = -M_planet / M_tot
    for j in eachindex(z)
        z[j] *= scale
    end
    return z
end
