# Generates the package logo.
#
# NOT PART OF THE DOCS BUILD. This is also the acceptance canary for the
# Makie extension: it exercises `orbit` construction with the M0+epoch phase
# group, `convert_arguments` quick-look plotting, `orbit_track_epochs`, and
# plain-Makie composition on top of them.
#
# Run from an environment with PlanetOrbits + CairoMakie, from the repo root:
#   julia --project=<env with CairoMakie> docs/logo.jl
#
# The orbital elements are the v0.11 logo's: v0.11's `τ` (fraction of a period past
# tref = 58849) is spelled as the equivalent mean anomaly at that epoch,
# `M0 = 2πτ`.

using PlanetOrbits
using PlanetOrbits: orbit
using CairoMakie

logocolors = Makie.Colors.JULIA_LOGO_COLORS

orbit1 = orbit(
    a=0.8, e=0.0, i=0.0, ω=0.0, Ω=0.0,
    M0=2π * 0.7, epoch=58849.0,
    plx=1000, M=1.0,
)
orbit2 = orbit(
    a=1.269, e=0.16, i=0.0, ω=120, Ω=0.0,
    M0=2π * 0.8, epoch=58849.0,
    plx=1000, M=1.0,
)
moon = orbit(
    a=0.274, e=0.0, i=0.0, ω=120, Ω=0.0,
    M0=2π * 0.0, epoch=58849.0,
    plx=1000, M=1.0,
)

t0 = 58849.0
sol1 = orbitsolve(orbit1, t0)
sol2 = orbitsolve(orbit2, t0)

fig = Figure(size=(300, 300), backgroundcolor=:transparent)
ax = Axis(fig[1, 1];
    limits=((-1600, 1300), (-1300, 1600)),
    backgroundcolor=:transparent,
    aspect=DataAspect())
hidedecorations!(ax)
hidespines!(ax)

# Star
scatter!(ax, [0], [0]; color=logocolors.blue, markersize=40,
    strokewidth=1, strokecolor="#222")

# Planet 1: quick-look convert_arguments path — `lines!(ax, sys)` of a
# two-body system traces the orbit track.
lines!(ax, orbit1; color=logocolors.green, linewidth=2.5)
scatter!(ax, [raoff(sol1)], [decoff(sol1)]; color=logocolors.green,
    markersize=26, strokewidth=1, strokecolor="#222")

# Planet 2
lines!(ax, orbit2; color=logocolors.red, linewidth=2.5)
x2, y2 = raoff(sol2), decoff(sol2)
scatter!(ax, [x2], [y2]; color=logocolors.red,
    markersize=18, strokewidth=1, strokecolor="#222")

# Moon of planet 2: the track drawn about planet 2's position.
ts_moon = orbit_track_epochs(moon; n=100)
traj_moon = orbitsolve(moon, ts_moon)
xs = raoff.(traj_moon) .+ x2
ys = decoff.(traj_moon) .+ y2
lines!(ax, xs, ys; color=logocolors.purple, linewidth=2.0)
solm = orbitsolve(moon, t0)
scatter!(ax, [raoff(solm) + x2], [decoff(solm) + y2]; color=logocolors.purple,
    markersize=12, strokewidth=1, strokecolor="#222")

save(joinpath(@__DIR__, "src", "assets", "logo.svg"), fig)
save(joinpath(@__DIR__, "src", "assets", "logo.png"), fig)
fig
