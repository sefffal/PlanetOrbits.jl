precompile(VisualElements, (NTuple{8,Float64},))
precompile(VisualElementsDeg, (NTuple{8,Float64},))
precompile(Core.kwfunc(VisualElements), (Vector{Any}, typeof(VisualElements), Int))
precompile(Core.kwfunc(VisualElementsDeg), (Vector{Any}, typeof(VisualElementsDeg), Int))

for f in (kep2cart, xyz, x, y, z)
    precompile(f, (VisualElements{Float64}, Float64))
    precompile(f, (VisualElements{Float64}, Int))
    precompile(f, (VisualElementsDeg{Float64}, Float64))
    precompile(f, (VisualElementsDeg{Float64}, Int))
end