precompile(KeplerianElements, (NTuple{8,Float64},))
precompile(KeplerianElementsDeg, (NTuple{8,Float64},))
precompile(Core.kwfunc(KeplerianElements), (Vector{Any}, typeof(KeplerianElements), Int))
precompile(Core.kwfunc(KeplerianElementsDeg), (Vector{Any}, typeof(KeplerianElementsDeg), Int))

for f in (xyz, x, y, z)
    precompile(f, (KeplerianElements{Float64}, Float64))
    precompile(f, (KeplerianElements{Float64}, Int))
    precompile(f, (KeplerianElementsDeg{Float64}, Float64))
    precompile(f, (KeplerianElementsDeg{Float64}, Int))
end