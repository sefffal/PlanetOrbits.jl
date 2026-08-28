# Matrix-based replica of `ahl21_step`, for high-precision references.
#
# `ahl21_step` builds MMatrix scratch, which StaticArrays refuses for a
# non-isbits eltype — so BigFloat and Dual{BigFloat} cannot go through it.
# This is the same map over plain Matrix, calling the *same* `_delxv_gamma` /
# `_comp_sum` for everything numerically substantive. Callers assert it
# reproduces `ahl21_step` bit-for-bit in Float64 before trusting it.
#
# Shared by perf/nbody-gradient-endtoend.jl and perf/nbody-enzyme.jl.

using StaticArrays
using PlanetOrbits: _delxv_gamma, _comp_sum, _dot3

function ref_step(x::Matrix{T}, v::Matrix{T}, xerr::Matrix{T}, verr::Matrix{T},
                  Gm::Vector{T}, h::T) where {T}
    N = length(Gm)
    drift!(hh) = for i in 1:N, k in 1:3
        x[k,i], xerr[k,i] = _comp_sum(x[k,i], xerr[k,i], hh * v[k,i])
    end
    function pair!(i, j, hh, drift_first)
        x0 = SVector{3,T}(x[1,i]-x[1,j], x[2,i]-x[2,j], x[3,i]-x[3,j])
        v0 = SVector{3,T}(v[1,i]-v[1,j], v[2,i]-v[2,j], v[3,i]-v[3,j])
        gm = Gm[i] + Gm[j]
        iszero(gm) && return
        delx, delv = _delxv_gamma(x0, v0, gm, hh, drift_first)
        gminv = inv(gm); mi = Gm[i]*gminv; mj = Gm[j]*gminv
        for k in 1:3
            x[k,i], xerr[k,i] = _comp_sum(x[k,i], xerr[k,i],  mj*delx[k])
            x[k,j], xerr[k,j] = _comp_sum(x[k,j], xerr[k,j], -mi*delx[k])
            v[k,i], verr[k,i] = _comp_sum(v[k,i], verr[k,i],  mj*delv[k])
            v[k,j], verr[k,j] = _comp_sum(v[k,j], verr[k,j], -mi*delv[k])
        end
    end
    function phi!(hh)
        a = zeros(T, 3, N)
        for i in 1:N-1, j in i+1:N
            rij = SVector{3,T}(x[1,i]-x[1,j], x[2,i]-x[2,j], x[3,i]-x[3,j])
            r2 = _dot3(rij); r3 = r2 * sqrt(r2)
            for k in 1:3
                fac = rij[k] / r3
                a[k,i] -= Gm[j] * fac
                a[k,j] += Gm[i] * fac
            end
        end
        coeff = hh^3 / 24
        for i in 1:N-1, j in i+1:N
            aij = SVector{3,T}(a[1,i]-a[1,j], a[2,i]-a[2,j], a[3,i]-a[3,j])
            rij = SVector{3,T}(x[1,i]-x[1,j], x[2,i]-x[2,j], x[3,i]-x[3,j])
            r2 = _dot3(rij); r1 = sqrt(r2)
            ardot = _dot3(aij, rij)
            fac1 = coeff / (r2 * r2 * r1)
            fac2 = 2 * (Gm[i] + Gm[j]) / r1 + 3 * ardot
            for k in 1:3
                fac = fac1 * (rij[k] * fac2 - r2 * aij[k])
                v[k,i], verr[k,i] = _comp_sum(v[k,i], verr[k,i],  Gm[j]*fac)
                v[k,j], verr[k,j] = _comp_sum(v[k,j], verr[k,j], -Gm[i]*fac)
            end
        end
    end
    h2 = h / 2
    drift!(h2)
    for i in 1:N-1, j in i+1:N; pair!(i, j, h2, true); end
    phi!(h)
    for i in N-1:-1:1, j in N:-1:i+1; pair!(i, j, h2, false); end
    drift!(h2)
    return nothing
end

