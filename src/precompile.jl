precompile(VisualOrbit, (NTuple{8,Float64},))
precompile(VisualOrbitDeg, (NTuple{8,Float64},))
precompile(Core.kwfunc(VisualOrbit), (Vector{Any}, typeof(VisualOrbit), Int))
precompile(Core.kwfunc(VisualOrbitDeg), (Vector{Any}, typeof(VisualOrbitDeg), Int))

for f in (kep2cart, xyz, x, y, z)
    precompile(f, (VisualOrbit{Float64}, Float64))
    precompile(f, (VisualOrbit{Float64}, Int))
    precompile(f, (VisualOrbitDeg{Float64}, Float64))
    precompile(f, (VisualOrbitDeg{Float64}, Int))
end