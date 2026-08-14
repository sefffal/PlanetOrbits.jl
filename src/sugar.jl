# ---------------------------------------------------
# Two-body convenience layer
# ---------------------------------------------------

"""
    orbit(; M, a, e=0, i=0, ω=0, Ω=0, tp=0, plx=…, ra=…, …)

Construct a trivial two-body `System`: a primary of mass `M` in solar masses,
named `:A`, with a massless secondary named `:b` on the given orbit. This
reproduces the classic PlanetOrbits v1 `KepOrbit`/`Visual`/`AbsoluteVisual`
semantics, where `M` is the total mass and one-argument observables
(`raoff(sol)`, …) give the secondary relative to the primary.

Frame keywords are passed through to `System` (none, `plx`, or the full
absolute frame set).

This v1-compatibility convenience is deliberately *not* exported — opt in
with `using PlanetOrbits: orbit`. For a secondary with real mass — and for
anything hierarchical — construct `System(Orbit(…))` explicitly.
"""
function orbit(; M, a=nothing, P=nothing, e=0.0, i=0.0, ω=0.0, Ω=0.0,
               tp=nothing, M0=nothing, θ=nothing, epoch=nothing, kwargs...)
    A = Body(mass=M, name=:A)
    b = Body(mass=zero(float(M)), name=:b)
    # The whole phase group forwards to `Orbit`; default tp=0 preserves the
    # v1 semantics when no phase parametrization is given.
    if tp === nothing && M0 === nothing && θ === nothing
        tp = 0.0
    end
    return System((A, b), (Orbit(b, about=A; a, P, e, i, ω, Ω, tp, M0, θ, epoch),); kwargs...)
end

"""
    ThieleInnes(; A, B, F, G, plx=nothing)

Convert Thiele-Innes constants to the size and orientation elements, as a
NamedTuple to splat into `Orbit`:

    Orbit(b, about=A; ThieleInnes(A=…, B=…, F=…, G=…, plx=24.5)...,
                      e=0.3, tp=59000.0)

Thiele-Innes replaces `a`, `i`, `ω` and `Ω` jointly; supply the shape and
phase elements as usual. This parametrization has no coordinate singularity
at `e → 0` or `i → 0`, which is why Gaia-only astrometric fits use it.

`A, B, F, G` are in AU, or in mas if `plx` [mas] is given.

!!! note "The node is ambiguous by ±180°"
    `(ω, Ω)` and `(ω+π, Ω+π)` give *identical* Thiele-Innes constants — every
    term picks up two sign changes — so the inverse is genuinely two-valued.
    The two solutions have the same sky-plane track and opposite line-of-sight
    motion, which is the familiar astrometric node ambiguity: astrometry alone
    cannot tell the ascending node from the descending one. Radial velocities
    break the tie.

    This function returns the branch with `Ω ∈ [0, π)`. If you have RV data
    preferring the other node, use `(ω+π, Ω+π)`.

The inverse is `thieleinnes(sys[, k]; plx=nothing)`.
"""
function ThieleInnes(; A, B, F, G, plx=nothing)
    A, B, F, G = promote(float(A), float(B), float(F), float(G))
    # Inverted directly from the forward definitions
    #   A =  cosΩcosω − sinΩsinω cosi,   B =  sinΩcosω + cosΩsinω cosi
    #   F = −cosΩsinω − sinΩcosω cosi,   G = −sinΩsinω + cosΩcosω cosi
    # which give the exact identities
    #   A+G = a(1+cosi)cos(ω+Ω)   B−F =  a(1+cosi)sin(ω+Ω)
    #   A−G = a(1−cosi)cos(ω−Ω)   B+F = −a(1−cosi)sin(ω−Ω)
    # so the two half-angle sums are plain atan2s and (a, i) come from the
    # two moduli. No quadrant fixups and no near-zero cosine division —
    # v1's ThieleInnesOrbit needed both and carried documented π errors for
    # Ω ≥ π and for ω+Ω > 2π.
    sp = hypot(A + G, B - F)              # a(1 + cos i)
    sm = hypot(A - G, B + F)              # a(1 − cos i)
    ωpΩ = atan(B - F, A + G)
    ωmΩ = atan(-(B + F), A - G)
    ω = rem2pi((ωpΩ + ωmΩ) / 2, RoundDown)
    Ω = rem2pi((ωpΩ - ωmΩ) / 2, RoundDown)
    # Both half-angle sums are unchanged by shifting ω and Ω together by π, so
    # pick the conventional branch Ω ∈ [0, π). See the note above: this is a
    # real degeneracy, not a numerical fixup.
    # `_primal`: which of the two equivalent (ω, Ω) representations comes back
    # is a branch, and a branch taken on a `Dual` can differ between the value
    # and its perturbation at exactly Ω = π.
    if _primal(Ω) >= π
        Ω -= π
        ω = rem2pi(ω - π, RoundDown)
    end
    α = (sp + sm) / 2                     # semi-major axis, input units
    i = acos(clamp((sp - sm) / (sp + sm), -one(α), one(α)))
    a = plx === nothing ? α : α / float(plx)
    return (; a, i, ω, Ω)
end

"""
    thieleinnes(sys, k=1; plx=nothing)

Thiele-Innes constants `(A, B, F, G)` of hierarchy row `k`, in AU — or in mas
if `plx` [mas] is given. Inverse of [`ThieleInnes`](@ref).
"""
function thieleinnes(sys::System, k::Integer=0; plx=nothing)
    row = k == 0 ? _only_row(sys) : sys.rows[k]
    α = plx === nothing ? row.a : row.a * float(plx)
    A = α * (row.cosΩ * row.cosω - row.sinΩ * row.sinω * row.cosi)
    B = α * (row.sinΩ * row.cosω + row.cosΩ * row.sinω * row.cosi)
    F = α * (-row.cosΩ * row.sinω - row.sinΩ * row.cosω * row.cosi)
    G = α * (-row.sinΩ * row.sinω + row.cosΩ * row.cosω * row.cosi)
    return (; A, B, F, G)
end

"""
    rvorbit(; M, msini=0, a=…|P=…, e=0, ω=0, tp=0)

Construct a two-body `System` in the radial-velocity-only convention: no
parallax (so angular observables are unavailable rather than silently wrong),
`i = π/2` and `Ω = 0`.

Radial velocities constrain only `m·sin i`, never `m` and `i` separately, so
`msini` is exactly what the data measure — under the `i = π/2` convention the
secondary's mass *is* its m·sin i. Reading it as a true mass is a lower bound.

Like [`orbit`](@ref), this is an opt-in convenience: `using PlanetOrbits: rvorbit`.
For a full 3D fit, build the `System` directly and let the astrometry
constrain `i`.
"""
function rvorbit(; M, msini=0.0, a=nothing, P=nothing, e=0.0, ω=0.0, tp=0.0)
    A = Body(mass=M, name=:A)
    b = Body(mass=msini, name=:b)
    return System((A, b), (Orbit(b, about=A; a, P, e, ω, tp, i=π / 2, Ω=0.0),))
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
_period(row::Row) = row.hyperbolic ? oftype(row.n, Inf) : 2π / row.n * year2day_julian

"""
    totalmass(sys)

Total mass of every body in the system [M⊙].
"""
totalmass(sys::System) = sum(sys.masses)

# Orbital-element accessors: the value for the system's only orbit, or for
# hierarchy row `k`.
for (fn, field, desc) in (
    (:eccentricity, :e, "Eccentricity of the system's only orbit, or of hierarchy row `k`. Greater than 1 for unbound orbits."),
    (:inclination, :i, "Inclination [rad] of the system's only orbit, or of hierarchy row `k`."),
    (:semimajoraxis, :a, "Semi-major axis [AU] of the system's only orbit, or of hierarchy row `k`. Negative for unbound (hyperbolic) orbits."),
    (:meanmotion, :n, "Mean motion [rad / julian year] of the system's only orbit, or of hierarchy row `k`."),
    (:periastron, :tp, "Epoch of periastron passage `tp` [MJD] of the system's only orbit, or of hierarchy row `k`."))
    @eval begin
        @doc """
            $($(string(fn)))(sys)
            $($(string(fn)))(sys, k)

        $($desc)
        """
        function $fn end
        $fn(sys::System) = _only_row(sys).$field
        $fn(sys::System, k::Integer) = sys.rows[k].$field
    end
end

"""
    distance(sys)

Distance to the system [pc]. Requires a parallax or absolute frame.
"""
distance(sys::System{NB,NR,T,NoFrame}) where {NB,NR,T} =
    error("this system has no parallax; construct it with plx=… to define a distance")
distance(sys::System) = distance(sys.frame)

export period, totalmass, eccentricity, inclination, semimajoraxis,
       meanmotion, periastron, distance, ThieleInnes, thieleinnes
