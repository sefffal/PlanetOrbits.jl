### A Pluto.jl notebook ###
# v0.14.2

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 1793f5fa-7f30-4a1d-baef-b06436e1fc71
using PlutoUI, Revise, DirectOrbits, Plots, StaticArrays; theme(:dao)

# ╔═╡ 872e19d4-c071-40bb-a091-22e48e85e2a6
md"""
# Orbit Playground

A Pluto notebook for visualizing orbits
"""

# ╔═╡ 60c8b78b-1b70-42d0-8fe6-7e40b0cdd4a2
M = 1.0;

# ╔═╡ 9da00dbb-c645-4bb9-a26c-5f148efb36cd
plx = 1000.;

# ╔═╡ f1ef0015-d671-450f-80cf-dc6651460998
τ = 0.0

# ╔═╡ 593a177f-bf0b-4b05-9408-745793ab2536
# begin 
# 	p1 = plot(elem, label="orbit", lw=2)
# 	scatter!([raoff(elem, t)], [decoff(elem, t)], label="planet", ms=10)
# 	scatter!([0], [0], label="star", marker=:star, ms=10)
	
# 	xlabel!("ΔRA - mas")
# 	ylabel!("ΔDEC - mas")
	
# 	ts = range(0, period(elem)*2, length=200)
# 	rv = radvel.(elem, ts)
# 	p2 = plot(ts, rv, label="orbit")
# 	scatter!([t], [radvel(elem, t)], label="planet", ms=10)
# 	xlabel!("t - days")
# 	ylabel!("Planet RV - km/s")
	
# 	layout = @layout [
# 		A{0.8h}
# 		B
# 	]
# 	plot(p1, p2; layout, size=(650,650))
# end

# ╔═╡ 994b3fb3-89d4-46e4-85de-041b9972ea91
@bind dt Slider(0:400, default=0)

# ╔═╡ c179bd9e-e392-4622-92a5-d64f442e2840
md"""  
a $(@bind a Slider(0.2:0.01:4, default=1))
i $(@bind i Slider(0:180.0, default=0))
e $( @bind e Slider(0:0.01:0.9))


t $(@bind t Slider(0:1:650, default=0))
Ω $(@bind Ω Slider(0:360.0, default=0))
ω $(@bind ω Slider(0:360.0, default=0))
"""

# ╔═╡ 3f38e8ca-2286-427f-8816-3b8b6cc78c74
elem = KeplerianElementsDeg(a,i,e,τ,M,ω,Ω,plx)

# ╔═╡ bfc82a87-4672-41fe-b136-ce15c84648d6
begin
	# points = [
	# 	0 1000
	# 	1000 0
	# 	0 -1000
	# 	-1000 0
	# ]
	# xs = points[:,1]
	# ys = points[:,2]
	ts = range(0, period(elem), length=20)
	xs = raoff.(elem, ts)
	ys = decoff.(elem, ts)
end

# ╔═╡ 24f48bae-f55a-4447-b5a9-c9637d03ef49
i

# ╔═╡ 7ab306ff-a7bc-4a7c-981d-9170a5398c02
ot = OrbitalTransformation(;i=deg2rad(i), e=Float64(e), M=Float64(M), ω=deg2rad(ω), Ω=deg2rad(Ω), plx=Float64(plx), platescale=1., dt=Float64(dt))

# ╔═╡ 9cd6b8c7-0608-4cc1-9037-342c436686c4
ot(SVector(0, 1000))

# ╔═╡ ba44719b-861a-4f19-b6ab-75717c37ee5e
begin
	mapped = [ot(SVector(x,y)) for (x,y) in zip(xs,ys)]
	newxs = [m[1] for m in mapped]
	newys = [m[2] for m in mapped]
end

# ╔═╡ d2c3bf24-fbdb-429d-96a2-b05ded00ae39
begin
	plot()
	scatter!(xs, ys)
	scatter!(newxs, newys, label="transformed")
	xlims!(-1200,1200)
	ylims!(-1200,1200)
end

# ╔═╡ Cell order:
# ╟─872e19d4-c071-40bb-a091-22e48e85e2a6
# ╠═1793f5fa-7f30-4a1d-baef-b06436e1fc71
# ╠═60c8b78b-1b70-42d0-8fe6-7e40b0cdd4a2
# ╠═9da00dbb-c645-4bb9-a26c-5f148efb36cd
# ╠═f1ef0015-d671-450f-80cf-dc6651460998
# ╠═593a177f-bf0b-4b05-9408-745793ab2536
# ╠═3f38e8ca-2286-427f-8816-3b8b6cc78c74
# ╠═24f48bae-f55a-4447-b5a9-c9637d03ef49
# ╠═7ab306ff-a7bc-4a7c-981d-9170a5398c02
# ╠═9cd6b8c7-0608-4cc1-9037-342c436686c4
# ╠═bfc82a87-4672-41fe-b136-ce15c84648d6
# ╠═994b3fb3-89d4-46e4-85de-041b9972ea91
# ╟─c179bd9e-e392-4622-92a5-d64f442e2840
# ╠═d2c3bf24-fbdb-429d-96a2-b05ded00ae39
# ╠═ba44719b-861a-4f19-b6ab-75717c37ee5e
