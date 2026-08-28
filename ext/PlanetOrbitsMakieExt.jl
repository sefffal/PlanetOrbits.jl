# ---------------------------------------------------
# PlanetOrbitsMakieExt
#
# Three things:
#
#   1. `MJDConversion` — a Makie dimension conversion so epoch axes carry
#      MJD Float64 data but show calendar-date ticks that recompute on every
#      zoom/pan (Makie's own DateTimeTicks locator does the work). Date and
#      DateTime values plot straight into the same axis.
#
#   2. `orbit_theme` — the house style (no grid lines; everything else is
#      deliberately left to Makie defaults so user themes compose).
#
#   3. Orbit-track plotting: `convert_arguments` so `lines(sys, :b, :A)`
#      works, and `orbitlines!` for phase-coloured tracks. Point spacing is
#      uniform eccentric anomaly via `orbit_track_epochs` (see plotting.jl).
#
# Everything here goes through the one public solve path (`orbitsolve`), so
# it works under either propagator.
# ---------------------------------------------------
module PlanetOrbitsMakieExt

using PlanetOrbits
using PlanetOrbits: System, Body, Row, NoFrame, _only_row, _setnames, _names,
                    _period, plotinfo
using Makie
using Dates

# ---------------------------------------------------
# MJD ↔ DateTime, pure arithmetic.
#
# Deliberately *not* PlanetOrbits.mjd/mjd2date: those go through AstroTime's
# TT epoch and are not mutual inverses at the ~minute level. Tick placement
# needs an exactly invertible pair; the ~1 min TT-vs-UTC distinction is far
# below anything an axis label resolves.
# ---------------------------------------------------
const _MJD_EPOCH_MS = Dates.value(DateTime(1858, 11, 17))
_mjd2dt(x::Real) = DateTime(Dates.UTM(round(Int64, _MJD_EPOCH_MS + x * 86_400_000)))
_dt2mjd(dt::Union{Date,DateTime}) = (Dates.value(DateTime(dt)) - _MJD_EPOCH_MS) / 86_400_000

struct MJDConv <: Makie.AbstractDimConversion end
PlanetOrbits.MJDConversion() = MJDConv()

_tomjd(x::Real) = Float64(x)
_tomjd(x::Union{Date,DateTime}) = _dt2mjd(x)

Makie.convert_dim_value(::MJDConv, value::Real) = Float64(value)
Makie.convert_dim_value(::MJDConv, value::Union{Date,DateTime}) = _dt2mjd(value)
Makie.convert_dim_value(::MJDConv, values::AbstractArray) = _tomjd.(values)
# The 4-arg form is what the plot pipeline calls on plot arguments.
Makie.convert_dim_value(::MJDConv, attr, values, previous_values) = _tomjd.(values)
Makie.convert_dim_value(::MJDConv, attr, value::Union{Real,Date,DateTime}, previous_values) =
    _tomjd(value)

# Mesh-shaped plots — `band!`, `poly!`, `hspan!` — reach the conversion after
# their arguments have already been assembled into vertices, so what arrives
# here is a vector of `Point2`s (plus a face list), not the x values a user
# passed. Those coordinates are already in MJD; the conversion has nothing to
# do but hand them back. Without this a `band!` on an epoch axis fails with
# `no method matching _tomjd(::Point{2, Float64})` — which is how the
# radial-velocity Gaussian-process band found it.
Makie.convert_dim_value(::MJDConv, values::AbstractArray{<:Makie.VecTypes}) = values
Makie.convert_dim_value(::MJDConv, attr, values::AbstractArray{<:Makie.VecTypes},
                        previous_values) = values
Makie.convert_dim_value(::MJDConv, value::Makie.VecTypes) = value
Makie.convert_dim_value(::MJDConv, attr, value::Makie.VecTypes, previous_values) = value

function Makie.get_ticks(::MJDConv, ticks, scale, formatter, vmin, vmax)
    if ticks isa Makie.Automatic && scale === identity
        # Delegate placement + neighbour-aware labelling to Makie's DateTime
        # machinery, then map the chosen datetimes back to MJD positions.
        dts, labels = Makie.get_ticks(Makie.DateTimeTicks(), scale, formatter,
                                      _mjd2dt(vmin), _mjd2dt(vmax))
        return _dt2mjd.(dts), labels
    end
    # Explicit user ticks (or a nonlinear scale): treat values as plain MJD.
    return Makie.get_ticks(ticks, scale, formatter, vmin, vmax)
end

"""
    add_mjd_axis!(gp, ax; position=:bottom, label="MJD", ticklabelscale=0.8) -> Axis

Companion numeric-MJD axis for an epoch axis `ax` built with
`MJDConversion` (which shows calendar dates). Creates a decoration-only
`Axis` in the same layout cell `gp`, x-linked to `ax`.

The companion is created *after* `ax` and therefore draws over it, so its
background is transparent — without that it would paint a white rectangle
across everything already plotted in the cell.
"""
function PlanetOrbits.add_mjd_axis!(gp, ax::Makie.Axis; position::Symbol=:bottom,
                                    label="MJD", ticklabelscale=0.8)
    ax2 = Makie.Axis(gp;
        xlabel=label,
        xaxisposition=position,
        backgroundcolor=:transparent,
        xgridvisible=false, ygridvisible=false,
        xticklabelsize=ticklabelscale * Makie.to_value(ax.xticklabelsize),
    )
    Makie.hideydecorations!(ax2)
    Makie.hidespines!(ax2)
    ax2.xzoomlock = true
    ax2.xrectzoom = false
    Makie.linkxaxes!(ax, ax2)
    return ax2
end

# ---------------------------------------------------
# Theme
# ---------------------------------------------------
function PlanetOrbits.orbit_theme()
    return Makie.Theme(
        Axis=(; xgridvisible=false, ygridvisible=false),
        Colorbar=(; ticksvisible=false),
    )
end

# ---------------------------------------------------
# Orbit tracks
# ---------------------------------------------------

const _NameLike = Union{Symbol,Body}

# The row a track belongs to: the one whose exterior is exactly this body.
# (EA spacing and phase colour are per *row*; a body's natural track is the
# orbit that places it.)
function _row_for(sys::System{NB,NR}, target::_NameLike) where {NB,NR}
    NR == 1 && return 1
    tname = target isa Symbol ? target : PlanetOrbits._name(target)
    names = _names(sys)
    for k in 1:NR
        ext = _setnames(names, sys.specs[k].ext)
        if ext == (tname,)
            return k
        end
    end
    error("no hierarchy row has $(tname) as its exterior; pass `row=k` explicitly " *
          "(rows: $(join(("$k: $(_setnames(names, s.ext)) about $(_setnames(names, s.int))" for (k, s) in enumerate(sys.specs)), "; ")))")
end

_isangular(sys::System) = !(sys.frame isa NoFrame)

function _trackpoints(sys::System, target, ref; row=nothing, n=150, tstart=nothing, solvekw=(;))
    k = row === nothing ? _row_for(sys, target) : row
    ts = orbit_track_epochs(sys, k; n, tstart)
    traj = orbitsolve(sys, ts; solvekw...)
    if _isangular(sys)
        xs = raoff.(traj, Ref(target), Ref(ref))
        ys = decoff.(traj, Ref(target), Ref(ref))
    else
        xs = posx.(traj, Ref(target), Ref(ref))
        ys = posy.(traj, Ref(target), Ref(ref))
    end
    phases = orbit_phase.(Ref(sys), k, ts)
    return xs, ys, phases
end

# `lines(sys, :b, :A)` / `scatter(sys, b, A)` quick looks.
function Makie.convert_arguments(::Makie.PointBased, sys::System,
                                 target::_NameLike, ref::_NameLike)
    xs, ys, _ = _trackpoints(sys, target, ref)
    return (Makie.Point2.(xs, ys),)
end

# Single-argument form for trivial two-body systems, matching the v0.11
# convert_single_argument convenience (`lines(orbit(...))`).
function Makie.convert_arguments(::Makie.PointBased, sys::System{NB,1}) where {NB}
    names = _names(sys)
    ext = _setnames(names, sys.specs[1].ext)
    int = _setnames(names, sys.specs[1].int)
    (length(ext) == 1 && length(int) == 1) ||
        error("pass an explicit target and reference: lines(sys, target, ref)")
    xs, ys, _ = _trackpoints(sys, ext[1], int[1])
    return (Makie.Point2.(xs, ys),)
end

"""
    orbitlines!(ax, sys, target, ref; row=nothing, n=150, tstart=nothing,
                colormap=..., colorbyphase=true, solvekw=(;), kwargs...)

Phase-coloured closed orbit track of `target` relative to `ref`. Sets the
axis orientation for sky plots (RA increasing left) via `plotinfo`.
"""
function PlanetOrbits.orbitlines!(ax::Makie.Axis, sys::System, target::_NameLike,
                                  ref::_NameLike;
                                  row=nothing, n=150, tstart=nothing,
                                  colormap=Makie.cgrad([Makie.wong_colors()[1], "#DDDDDD"]),
                                  colorbyphase=true,
                                  solvekw=(;),
                                  kwargs...)
    xs, ys, phases = _trackpoints(sys, target, ref; row, n, tstart, solvekw)
    if _isangular(sys) && !ax.xreversed[]
        ax.xreversed = true
    end
    if colorbyphase
        return Makie.lines!(ax, xs, ys; color=phases, colormap, colorrange=(0, 2π), kwargs...)
    else
        return Makie.lines!(ax, xs, ys; kwargs...)
    end
end

end # module
