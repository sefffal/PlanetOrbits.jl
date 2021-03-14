precompile(KeplarianElements, (NTuple{8,Float64},))
precompile(Core.kwfunc(KeplarianElements), (Vector{Any}, typeof(KeplarianElements), Int))

for f in (xyz, x, y, z)
    precompile(f, (KeplarianElements{Float64}, Float64))
    precompile(f, (KeplarianElements{Float64}, Int))
end