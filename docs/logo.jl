using DirectOrbits
using Plots
using Colors

logocolors = Colors.JULIA_LOGO_COLORS

orbit1 = KeplerianElementsDeg(
    a = 0.8,
    i = 0.0,
    e = 0.0,
    ω = 0.0,
    Ω = 0.0,
    τ = 0.7,
    plx=1000,
    μ = 1.0,
)

orbit2 = KeplerianElementsDeg(
    a = 1.4,
    i = 0.0,
    e = 0.16,
    ω = 120,
    Ω = 0.0,
    τ = 0.2,
    plx=1000,
    μ = 1.0,
)

p = plot(xlims=(-1600,1300), ylims=(-1300,1600), size=(300,300), framestyle=:none, legend=nothing, margin=0Plots.mm, background=:transparent)
scatter!([0], [0], color=logocolors.blue, markersize=20, markerstrokewidth=1, markerstrokecolor="#222")

plot!(orbit1, color=logocolors.green, linewidth=2.5)
scatter!([raoff(orbit1, 0)], [decoff(orbit1, 0)], color=logocolors.green, markersize=13, markerstrokewidth=1, markerstrokecolor="#222")

plot!(orbit2, color=logocolors.red, linewidth=2.5)
x = raoff(orbit2, 0)
y = decoff(orbit2, 0)
scatter!([x],[y], color=logocolors.red, markersize=9, markerstrokewidth=1, markerstrokecolor="#222")


moon = KeplerianElementsDeg(
    a = 0.2,
    i = 0,
    e = 0.0,
    ω = 120,
    Ω = 0.0,
    τ = 0.0,
    plx=1000,
    μ = 1.0,
)
νs = range(0, 2π, length=100)
xs = getproperty.(DirectOrbits.kep2cart_ν.(moon, νs), :x) .+ x
ys = getproperty.(DirectOrbits.kep2cart_ν.(moon, νs), :y) .+ y
plot!(xs,ys, color=logocolors.purple, linewidth=2.5)
i = 2
scatter!(xs[i:i],ys[i:i], color=logocolors.purple, markersize=6, markerstrokewidth=1, markerstrokecolor="#222")
savefig("docs/src/assets/logo.svg")
savefig("docs/src/assets/logo.png")
p