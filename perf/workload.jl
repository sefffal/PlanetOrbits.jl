# Shared definition of the Gaia-DR4-shaped benchmark workload: 150 epochs of
# absolute astrometry, 1–3 companions, and one likelihood-eval's worth of work
# (rebuild the system from θ, solve every epoch, read out the star's absolute
# sky path + proper motion + RV).
#
# Included by `gaia-workload.jl` (the §11 headline gate) and `dual-passes.jl`
# (the §10.4.1 per-pass gradient decomposition) so both measure the same model.

using PlanetOrbits
using PlanetOrbits: Body, Orbit, System

const NEP = 150
const EPOCHS = collect(range(57000.0, 62000.0, length=NEP))

# θ layout: [M_star, (m, a, e, i, ω, Ω, tp) per companion..., plx, ra, dec,
# pmra, pmdec, rv, ref_epoch]
function example_theta(NP::Int)
    θ = [1.1]
    as = (2.5, 8.0, 17.0)
    for p in 1:NP
        append!(θ, [3mjup + p * mjup, as[p], 0.1 + 0.1p, 0.5, 1.1, 2.2, 58849.0 - 300p])
    end
    append!(θ, [24.5, 45.0, -30.0, 100.0, -50.0, 25e3, 57388.5])
    return θ
end

@inline _frame_kw(θ, NP) = (
    plx=θ[2+7NP], ra=θ[3+7NP], dec=θ[4+7NP], pmra=θ[5+7NP],
    pmdec=θ[6+7NP], rv=θ[7+7NP], ref_epoch=θ[8+7NP])

@inline _planet(θ, p) = (mass=θ[7p-5], a=θ[7p-4], e=θ[7p-3], i=θ[7p-2],
                         ω=θ[7p-1], Ω=θ[7p], tp=θ[7p+1])

function build(θ, ::Val{1})
    A = Body(mass=θ[1], name=:A)
    b = Body(mass=_planet(θ, 1).mass, name=:b)
    q = _planet(θ, 1)
    System((A, b), (Orbit(b, about=A; a=q.a, e=q.e, i=q.i, ω=q.ω, Ω=q.Ω, tp=q.tp),);
        _frame_kw(θ, 1)...)
end
function build(θ, ::Val{2})
    A = Body(mass=θ[1], name=:A)
    q1 = _planet(θ, 1); q2 = _planet(θ, 2)
    b = Body(mass=q1.mass, name=:b); c = Body(mass=q2.mass, name=:c)
    System((A, b, c), (
        Orbit(b, about=A;      a=q1.a, e=q1.e, i=q1.i, ω=q1.ω, Ω=q1.Ω, tp=q1.tp),
        Orbit(c, about=(A, b); a=q2.a, e=q2.e, i=q2.i, ω=q2.ω, Ω=q2.Ω, tp=q2.tp));
        _frame_kw(θ, 2)...)
end
function build(θ, ::Val{3})
    A = Body(mass=θ[1], name=:A)
    q1 = _planet(θ, 1); q2 = _planet(θ, 2); q3 = _planet(θ, 3)
    b = Body(mass=q1.mass, name=:b); c = Body(mass=q2.mass, name=:c); d = Body(mass=q3.mass, name=:d)
    System((A, b, c, d), (
        Orbit(b, about=A;         a=q1.a, e=q1.e, i=q1.i, ω=q1.ω, Ω=q1.Ω, tp=q1.tp),
        Orbit(c, about=(A, b);    a=q2.a, e=q2.e, i=q2.i, ω=q2.ω, Ω=q2.Ω, tp=q2.tp),
        Orbit(d, about=(A, b, c); a=q3.a, e=q3.e, i=q3.i, ω=q3.ω, Ω=q3.Ω, tp=q3.tp));
        _frame_kw(θ, 3)...)
end

# The per-sample readout: what a Gaia IAD + HGCA + RV joint fit reads out.
function readout(traj, sys, ::Type{T}) where {T}
    refs = bodies(sys)
    bary = barycentre(sys)
    star = refs.A
    acc = zero(T)
    for sol in traj
        acc += raoff(sol, star, bary) + decoff(sol, star, bary) +
               pmra(sol, star, bary) + frame_pmra(sol) +
               pmdec(sol, star, bary) + frame_pmdec(sol) +
               radvel(sol, star, bary) + frame_rv(sol)
    end
    return acc
end

function evalonce(θ, np::Val{NP}, method, traj=nothing) where {NP}
    sys = build(θ, np)
    if traj === nothing
        traj = Trajectory{eltype(θ)}(sys, EPOCHS)
    end
    orbitsolve!(traj, sys; method)
    return readout(traj, sys, eltype(θ))
end
