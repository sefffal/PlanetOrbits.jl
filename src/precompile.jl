precompile(KeplarianElements, (NTuple{8,Float64},))
precompile(KeplarianElementsDeg, (NTuple{8,Float64},))
precompile(Core.kwfunc(KeplarianElements), (Vector{Any}, typeof(KeplarianElements), Int))
precompile(Core.kwfunc(KeplarianElementsDeg), (Vector{Any}, typeof(KeplarianElementsDeg), Int))

for f in (xyz, x, y, z)
    precompile(f, (KeplarianElements{Float64}, Float64))
    precompile(f, (KeplarianElements{Float64}, Int))
    precompile(f, (KeplarianElementsDeg{Float64}, Float64))
    precompile(f, (KeplarianElementsDeg{Float64}, Int))
end