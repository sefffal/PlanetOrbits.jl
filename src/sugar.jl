# ---------------------------------------------------
# Two-body convenience layer
# ---------------------------------------------------

"""
    orbit(; M, a, e=0, i=0, ω=0, Ω=0, tp=0, plx=…, ra=…, …)

Construct a trivial two-body `System`: a primary of mass `M` [M⊙] (named
`:A`) with a massless secondary (named `:b`) on the given orbit. This
reproduces the classic PlanetOrbits v1 `KepOrbit`/`Visual`/`AbsoluteVisual`
semantics, where `M` is the total mass and one-argument observables
(`raoff(sol)`, …) give the secondary relative to the primary.

Frame keywords are passed through to `System` (none, `plx`, or the full
absolute frame set).

This v1-compatibility convenience is deliberately *not* exported — opt in
with `using PlanetOrbits: orbit`. For a secondary with real mass — and for
anything hierarchical — construct `System(Orbit(…))` explicitly.
"""
function orbit(; M, a, e=0.0, i=0.0, ω=0.0, Ω=0.0, tp=0.0, kwargs...)
    A = Body(mass=M, name=:A)
    b = Body(mass=zero(float(M)), name=:b)
    return System(Orbit(b, about=A; a, e, i, ω, Ω, tp); kwargs...)
end

# ---------------------------------------------------
# System / orbit property accessors
#
# Row-index forms work on any system (rows are in depth-first post-order,
# innermost binaries first). One-argument forms require a single-row system.
# ---------------------------------------------------

@inline function _only_row(sys::System{NB,NR}) where {NB,NR}
    NR == 1 || error("this system has $NR orbits: pass a row index, e.g. period(sys, 1)")
    return sys.rows[1]
end

"""
    period(sys)
    period(sys, k)

Period [days] of the system's only orbit, or of hierarchy row `k`.
"""
period(sys::System) = _period(_only_row(sys))
period(sys::System, k::Integer) = _period(sys.rows[k])
_period(row::Row) = 2π / row.n * year2day_julian

totalmass(sys::System) = sum(sys.masses)

# eccentricity(sys[, k]), inclination(sys[, k]), semimajoraxis(sys[, k]),
# meanmotion(sys[, k]), periastron(sys[, k]):
# orbital-element accessors of the only orbit, or of hierarchy row `k`.
# meanmotion is rad/julian year; periastron is the periastron epoch tp [MJD].
for (fn, field) in ((:eccentricity, :e), (:inclination, :i),
                    (:semimajoraxis, :a), (:meanmotion, :n), (:periastron, :tp))
    @eval $fn(sys::System) = _only_row(sys).$field
    @eval $fn(sys::System, k::Integer) = sys.rows[k].$field
end

"""
    distance(sys)

Distance to the system [pc]. Requires a parallax or absolute frame.
"""
distance(sys::System{NB,NR,T,NoFrame}) where {NB,NR,T} =
    error("this system has no parallax; construct it with plx=… to define a distance")
distance(sys::System) = distance(sys.frame)

export period, totalmass, eccentricity, inclination, semimajoraxis,
       meanmotion, periastron, distance
