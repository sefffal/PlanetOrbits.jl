# ---------------------------------------------------
# Plotting support that needs no plotting package
#
# Two things live here so that PlanetOrbits and every downstream package
# (Octofitter's plot layer, user scripts) share one implementation:
#
#   1. `plotinfo` — the resolver table mapping each observable function to
#      its axis label, unit, and axis conventions (RA axes flip, position
#      angles wrap). v1 carried this as a NamedTuple inside a Plots.jl
#      recipe; it is package-agnostic metadata about the observables, so it
#      belongs beside them.
#
#   2. The epoch-grid heuristics — where to sample an orbit so a plotted
#      track looks smooth with few points. Uniform *eccentric anomaly* is
#      the right spacing for a closed track (dense at periastron, where the
#      curvature is); Kepler's equation is closed-form in that direction
#      (t from E needs no solver), so the grid is generated as plain epochs
#      and solved through the one normal `orbitsolve` path, under either
#      propagator.
#
# The Makie extension builds on these; it adds nothing that other backends
# could not also use.
# ---------------------------------------------------

"""
    plotinfo(f) -> (; label, unit, flip, wrap)

Axis metadata for an observable function `f` (e.g. `raoff`, `radvel`):

  - `label` — human-readable quantity name.
  - `unit`  — unit string, or `""` for dimensionless.
  - `flip`  — whether the axis should increase leftward/downward (RA-like
    axes: east is to the left on the sky).
  - `wrap`  — the period the quantity wraps with (`2π` for position angle),
    or `nothing`.

`plotlabel(f)` renders these as `"label [unit]"`.
"""
plotinfo(@nospecialize f) = (; label=string(f), unit="", flip=false, wrap=nothing)

plotinfo(::typeof(posx)) = (; label="Δx (east)", unit="au", flip=true, wrap=nothing)
plotinfo(::typeof(posy)) = (; label="Δy (north)", unit="au", flip=false, wrap=nothing)
plotinfo(::typeof(posz)) = (; label="Δz (line of sight)", unit="au", flip=false, wrap=nothing)
plotinfo(::typeof(velx)) = (; label="Δvx (east)", unit="au/yr", flip=true, wrap=nothing)
plotinfo(::typeof(vely)) = (; label="Δvy (north)", unit="au/yr", flip=false, wrap=nothing)
plotinfo(::typeof(velz)) = (; label="Δvz (line of sight)", unit="au/yr", flip=false, wrap=nothing)
plotinfo(::typeof(radvel)) = (; label="radial velocity", unit="m/s", flip=false, wrap=nothing)
plotinfo(::typeof(raoff)) = (; label="Δα*", unit="mas", flip=true, wrap=nothing)
plotinfo(::typeof(decoff)) = (; label="Δδ", unit="mas", flip=false, wrap=nothing)
plotinfo(::typeof(pmra)) = (; label="Δμα*", unit="mas/yr", flip=true, wrap=nothing)
plotinfo(::typeof(pmdec)) = (; label="Δμδ", unit="mas/yr", flip=false, wrap=nothing)
plotinfo(::typeof(projectedseparation)) = (; label="separation", unit="mas", flip=false, wrap=nothing)
plotinfo(::typeof(posangle)) = (; label="position angle", unit="rad", flip=false, wrap=2π)

"""
    plotlabel(f) -> String

Axis label `"label [unit]"` for an observable function, from [`plotinfo`](@ref).
"""
function plotlabel(@nospecialize f)
    info = plotinfo(f)
    return isempty(info.unit) ? String(info.label) : "$(info.label) [$(info.unit)]"
end

"""
    paraminfo(name::Symbol) -> (; label, unit, angle) | nothing

Display metadata for a *parameter* by its conventional name — the orbital
element keywords `Orbit` accepts, plus frame variables and common fit
parameters. `angle` marks radian-valued parameters conventionally displayed
in degrees. Corner plots and summary tables share this table instead of
each keeping a private label dictionary; `nothing` means "not a name this
table knows".
"""
function paraminfo(name::Symbol)
    haskey(_PARAMINFO, name) || return nothing
    return _PARAMINFO[name]
end

const _PARAMINFO = Dict{Symbol,NamedTuple{(:label, :unit, :angle),Tuple{String,String,Bool}}}(
    :a      => (label="semi-major axis", unit="au", angle=false),
    :P      => (label="period", unit="days", angle=false),
    :e      => (label="eccentricity", unit="", angle=false),
    :ω      => (label="argument of periapsis", unit="°", angle=true),
    :secosω => (label="√e cos ω", unit="", angle=false),
    :sesinω => (label="√e sin ω", unit="", angle=false),
    :ecosω  => (label="e cos ω", unit="", angle=false),
    :esinω  => (label="e sin ω", unit="", angle=false),
    :i      => (label="inclination", unit="°", angle=true),
    :Ω      => (label="position angle of\nascending node", unit="°", angle=true),
    :tp     => (label="epoch of periastron\npassage", unit="mjd", angle=false),
    :M0     => (label="mean anomaly\nat ref. epoch", unit="°", angle=true),
    :θ      => (label="position angle\nat ref. epoch", unit="°", angle=true),
    # v2 masses are uniformly M⊙ (they multiply GM_sun); v1's Mⱼᵤₚ planet
    # convention does not carry over.
    :mass   => (label="mass", unit="M⊙", angle=false),
    :M      => (label="total mass", unit="M⊙", angle=false),
    :plx    => (label="parallax", unit="mas", angle=false),
    :ra     => (label="RA", unit="°", angle=false),
    :dec    => (label="Dec", unit="°", angle=false),
    :pmra   => (label="μα*", unit="mas/yr", angle=false),
    :pmdec  => (label="μδ", unit="mas/yr", angle=false),
    :rv     => (label="RV", unit="m/s", angle=false),
)

# ---------------------------------------------------
# Epoch grids
# ---------------------------------------------------

"""
    orbit_track_epochs(sys, k=only; n=150, tstart=nothing) -> Vector{Float64}

Epochs [MJD] tracing one full cycle of hierarchy row `k`, spaced uniformly
in **eccentric anomaly** — dense at periastron, sparse at apastron — so the
plotted track is smooth with few points. The first and last epoch describe
the same orbital phase one period apart, closing the track.

The cycle drawn is the one containing `tstart` (default: the cycle starting
at the row's periastron epoch). Under `KeplerianApprox` every cycle is
identical; under time-dependent effects (absolute frames, N-body) the choice
of cycle matters, so pass `tstart` near your data.

Errors for hyperbolic rows — an unbound orbit has no period; sample a time
span instead.
"""
orbit_track_epochs(sys::System; kwargs...) = _track_epochs(_only_row(sys); kwargs...)
orbit_track_epochs(sys::System, k::Integer; kwargs...) = _track_epochs(sys.rows[k]; kwargs...)

function _track_epochs(row::Row; n::Integer=150, tstart=nothing)
    row.hyperbolic && error(
        "this orbit is unbound (e ≥ 1): it has no period, so a closed track is undefined. " *
        "Sample a time span instead, e.g. range(t0, t1, length=n).")
    P = _period(row)
    e = row.e
    # Kepler's equation forward: M = E − e sin E, t = tp + M/n. No solver.
    Es = range(0.0, 2π, length=Int(n))
    ts = [row.tp + (E - e * sin(E)) / (2π) * P for E in Es]
    if tstart !== nothing
        # Shift by a whole number of periods so the track covers the cycle
        # containing tstart.
        ts .+= floor((tstart - row.tp) / P) * P
    end
    return ts
end

"""
    plot_epochs(sys, tmin, tmax; points_per_period=30, min_points=200, max_points=1000)

Epochs [MJD] for plotting model curves over the span `[tmin, tmax]`,
allocated **per hierarchy row**: each bound row contributes points at
`period/points_per_period` spacing over the whole span, so a short-period
inner orbit stays smooth without inflating the grid for a long-period outer
one. A uniform base grid of `min_points` guarantees coverage; the union is
clamped to `max_points` by thinning the densest contributions first (the
uniform base and the data span are never dropped).

This replaces sampling the whole span at the *minimum* period over all
orbits, which over-resolves every wide orbit in a system with a broad period
spread.
"""
function plot_epochs(sys::System{NB,NR}, tmin::Real, tmax::Real;
                     points_per_period::Integer=30,
                     min_points::Integer=200,
                     max_points::Integer=1000) where {NB,NR}
    tmin < tmax || error("need tmin < tmax; got $tmin, $tmax")
    span = float(tmax - tmin)
    ts = collect(range(tmin, tmax, length=min_points))
    for k in 1:NR
        row = sys.rows[k]
        row.hyperbolic && continue
        P = _period(row)
        step = P / points_per_period
        step >= span && continue                     # base grid already covers it
        npts = min(ceil(Int, span / step), 100_000)  # guard degenerate tiny periods
        append!(ts, range(tmin, tmax, length=npts))
    end
    sort!(unique!(ts))
    if length(ts) > max_points
        # Thin uniformly by index; keeps relative density while bounding cost.
        idx = round.(Int, range(1, length(ts), length=max_points))
        ts = ts[unique(idx)]
    end
    return ts
end

"""
    orbit_phase(row_or_sys[, k], t) -> Float64

Mean-anomaly phase in `[0, 2π)` of hierarchy row `k` at epoch `t` [MJD],
computed from the row's elements (`2π·(t − tp)/P`, wrapped). This is what
plots colour orbit tracks by; it is exact under `KeplerianApprox` and the
osculating phase under `AHL21`. Returns `NaN` for hyperbolic rows.
"""
orbit_phase(sys::System, t::Real) = _orbit_phase(_only_row(sys), t)
orbit_phase(sys::System, k::Integer, t::Real) = _orbit_phase(sys.rows[k], t)
function _orbit_phase(row::Row, t::Real)
    row.hyperbolic && return NaN
    return rem2pi(2π * (t - row.tp) / _period(row), RoundDown)
end

export plotinfo, plotlabel, paraminfo, orbit_track_epochs, plot_epochs, orbit_phase

# ---------------------------------------------------
# Makie-extension function stubs
#
# Defined here so the names exist (and are documentable) without Makie; the
# extension adds the working methods.
# ---------------------------------------------------

"""
    MJDConversion()

A Makie dimension conversion for epoch axes: data are Modified Julian Days
(plain `Float64`), tick labels are calendar dates that adapt as you zoom.
`Date`/`DateTime` values plot into the same axis and convert automatically
(for single-argument plots like `vlines!`, which bypass Makie's dim
conversions, use `mjd`):

    ax = Axis(fig[1,1]; dim1_conversion=PlanetOrbits.MJDConversion())
    lines!(ax, epochs_mjd, rvs)
    scatter!(ax, [DateTime("2035-06-01")], [0.0])
    vlines!(ax, mjd("2035-06-01"))

Requires Makie to be loaded (the method is added by the PlanetOrbitsMakie
extension).
"""
function MJDConversion end

"""
    orbit_theme()

A Makie theme carrying the PlanetOrbits/Octofitter plot look: no axis grid
lines, and no colorbar tick marks. Deliberately minimal — everything else is
left at Makie's defaults so that a user theme composes with it rather than
being overridden. Requires Makie to be loaded.

    with_theme(orbit_theme()) do
        octoplot(model, chain)
    end
"""
function orbit_theme end

"""
    orbitlines!(ax, sys, target, ref; n=150, kwargs...)

Draw the closed orbit track(s) of `target` relative to `ref` on `ax`,
sampled uniformly in eccentric anomaly and coloured by orbital phase.
Requires Makie to be loaded.
"""
function orbitlines! end

"""
    add_mjd_axis!(axisparent, ax) -> Axis

Add a companion x-axis displaying raw MJD values below/above an epoch axis
built with [`MJDConversion`](@ref). Requires Makie to be loaded.
"""
function add_mjd_axis! end

export MJDConversion, orbit_theme, orbitlines!, add_mjd_axis!
