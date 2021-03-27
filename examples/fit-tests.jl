using Optim
using Plots
using DirectOrbits
theme(:dao)
##


points = [
    0.0      500.   0.0     0.0
    # year2days/4    -700.0    -800.0  0.0
    DirectOrbits.year2days/2    -1500.0     0.0     0.0
    DirectOrbits.year2days/2    -1500.0     0.0     0.0
]

function fit(
    a,
    i,
    e,
    τ,
    μ,
    ω,
    Ω,
    plx,
)
    el = KeplerianElements(a, i, e, τ, μ, ω, Ω, plx)
    ssq = 0.0
    for i in 1:size(points,1)
        point = points[1,:]
        out = kep2cart(el, point[1])
        ssq += sum((out .- point[2:end]).^2)
    end
    return ssq
end

# initial = [1.0, 0.00, 0.3, 365/4, 1.0, 0.01, 0.01, 1000.0]
initial =       [0.9, 0, 0.3, 0.0, 1.0, deg2rad(90), 0.00, 1000.0]
initial_fixed = [0.9, 0, 0.3, 0.0, deg2rad(90), 0.00]


##
fitvec(args) = fit(args[1:4]..., 1.0, args[5:6]..., 1000.0)
fitvec(initial)
results = optimize(fitvec, initial_fixed, NelderMead(), Optim.Options(iterations =1000))


##
out = Optim.minimizer(results)
bestfit = KeplerianElements(out[1:4]..., 1.0, out[5:6]..., 1000.0)
display(bestfit)
plot(xlims=(-1700,1700), ylims=(-1700,1700),legend=:topleft)


plot!(KeplerianElements(initial...), label="Initial")
scatter!(points[:,2], points[:,3], ms=3, color=:black, label="Data", markerz=1:3)

plot!(bestfit, label="Best fit")

# scatter!(raoff.(bestfit, points[:,1]), decoff.(bestfit, points[:,1]), markerz=1:3) 

##

using Zygote
##
fitvec(args) = fit(args...)
fitvec′(out, args) = out .= Zygote.gradient(fit, args...)
od = OnceDifferentiable(fitvec, fitvec′, initial)

# fitvec′′(out, args) = out .= Zygote.hessian(fit, args...)
# od = TwiceDifferentiable(fitvec, fitvec′, fitvec′′, initial)
results = optimize(od, initial, GradientDescent())
##
bestfit = KeplarianElements(Optim.minimizer(results)...)
plot(bestfit)
scatter!(points[:,2], points[:,3], ms=3, color=:black)


