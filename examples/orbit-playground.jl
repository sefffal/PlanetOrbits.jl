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
using Revise, PlutoUI, DirectOrbits, Plots; theme(:dao)

# ╔═╡ 872e19d4-c071-40bb-a091-22e48e85e2a6
md"""
# Orbit Playground

A Pluto notebook for visualizing orbits
"""

# ╔═╡ ce86c343-b6e4-4833-935f-ee0392d1ee89
md"""
*Specify the Keplerian elements to visualize using the text fields and sliders below.*
"""

# ╔═╡ 79f416a6-7c52-4d57-8bdd-1f84dc04d7e8
md"Gravitational Parameter ($M_⊙$)"

# ╔═╡ 60c8b78b-1b70-42d0-8fe6-7e40b0cdd4a2
μ = 1.0;

# ╔═╡ 124025b9-f7de-4395-a5d3-1dc6ed4a67f7
md"Paralax Distance (mas)"

# ╔═╡ 9da00dbb-c645-4bb9-a26c-5f148efb36cd
plx = 45.;

# ╔═╡ 800b91ae-d6c6-479f-afd7-d07d7b207cd3
md"Epoch of periastron passage [0,1]"

# ╔═╡ f1ef0015-d671-450f-80cf-dc6651460998
τ = 0.;

# ╔═╡ c179bd9e-e392-4622-92a5-d64f442e2840
md"""
a $(@bind a Slider(1:0.01:50, default=15))
i $(@bind i Slider(0:180.0, default=0))
e $( @bind e Slider(0:0.01:0.9, default=0))


t $(@bind t_frac Slider(0:0.01:1, default=0))
Ω $(@bind Ω Slider(0:360.0, default=0))
ω $(@bind ω Slider(0:360.0, default=0))
"""

# ╔═╡ 6addb17b-1a15-4e66-98d7-03e287664c34
md"""
Plot apparent speed $(@bind show_app_speed CheckBox())

(*The chance of seeing a planet at any given point of its orbit with direct imaging, all things being equal, is proportional to it's apparent / projected speed across the sky*)
"""

# ╔═╡ 3f38e8ca-2286-427f-8816-3b8b6cc78c74
elem = KeplerianElementsDeg(a,i,e,τ,μ,ω,Ω,plx)

# ╔═╡ 596e2d59-203c-4e69-985b-f8a82624ef6c
md"Time range in MJD (modified Juian days). Increase `length` to increase the resolution of the plots."

# ╔═╡ 465e92c8-1004-47d6-ac4e-69172afad2b0
ts = 58849 .+ range(0, 2period(elem), length=200);

# ╔═╡ 6f7a36ed-2dd2-4d93-a4cb-7e8c46fbdf22
t = t_frac * (last(ts)-first(ts)) + first(ts)

# ╔═╡ e8d619bc-e37d-437a-96fb-0995aed2f823
begin
	posn = kep2cart.(elem, ts)
	ra = [p[1] for p in posn]
	dec = [p[2] for p in posn]
	rv = [p[4] for p in posn]
end;

# ╔═╡ 9dd26db3-e443-46f3-8e18-21eb37b4d5b6
begin
	dradt = diff(ra)./step(ts).*DirectOrbits.year2days
	ddecdt = diff(dec)./step(ts).*DirectOrbits.year2days
	app_speed = sqrt.(dradt.^2 .+ ddecdt.^2)
end;

# ╔═╡ 593a177f-bf0b-4b05-9408-745793ab2536
begin 
	p1 = plot(xflip=true, aspectratio=1, legend=:topleft)
	
	# Option for plotting the RV as the orbit line colour
	# plot!(ra, dec, linez=rv, colorbar=:none, label="orbit", color=:diverging_bkr_55_10_c35_n256)
	
	if show_app_speed
		clims=(NaN, NaN)
		if abs(-(extrema(app_speed)...)) < 1e-2
			clims=(-1,1) .+ first(app_speed)
		end
		plot!(ra, dec; linez=app_speed, color=:plasma, label="orbit", colorbar_title="mas/year", lw=8, clims)
	
	else
		plot!(elem, label="orbit", color=1)
	end
	
	
	scatter!([raoff(elem, t)], [decoff(elem, t)], label="planet", ms=10, color=2)
	scatter!([0], [0], label="star", marker=:star, color=3, ms=10)
	
	xlabel!("ΔRA - mas")
	ylabel!("ΔDEC - mas")
	
	p2 = plot(ts, rv, legend=:none)
	scatter!([t], [radvel(elem, t)], ms=10)
	xlabel!("t - mjd")
	ylabel!("\$\\mathrm{RV_{planet} - km/s}\$")
	
	layout = @layout [
		A{0.8h}
		B
	]
	plot(p1, p2; layout, size=(650,650))
end

# ╔═╡ 36bc055e-3b5b-41a0-863a-53b78a6328d9
md"""
### CSV Export
Click below to download a CSV with the ΔRA, ΔDEC, and RV values
for t = $(round(Int, first(t))) to $(round(Int, first(t))) (mjd)
"""

# ╔═╡ 30c618b6-bd58-4335-bc55-23c16317011d
let
	# Create a super simple CSV export
	csv = "dRA (mas), dDEC (mas), RV (km/s), App. Speed (mas/year)\n"*join((
		join(p, ", ")
		for p in zip(ra,dec,rv,app_speed)
	), "\n")
	DownloadButton(csv, "orbit-a-$a-i-$i-e-$e-ω-$ω-Ω-$Ω-τ-$τ-μ-$μ.csv")
end

# ╔═╡ Cell order:
# ╟─872e19d4-c071-40bb-a091-22e48e85e2a6
# ╠═1793f5fa-7f30-4a1d-baef-b06436e1fc71
# ╟─ce86c343-b6e4-4833-935f-ee0392d1ee89
# ╟─79f416a6-7c52-4d57-8bdd-1f84dc04d7e8
# ╠═60c8b78b-1b70-42d0-8fe6-7e40b0cdd4a2
# ╟─124025b9-f7de-4395-a5d3-1dc6ed4a67f7
# ╠═9da00dbb-c645-4bb9-a26c-5f148efb36cd
# ╟─800b91ae-d6c6-479f-afd7-d07d7b207cd3
# ╠═f1ef0015-d671-450f-80cf-dc6651460998
# ╟─c179bd9e-e392-4622-92a5-d64f442e2840
# ╟─6f7a36ed-2dd2-4d93-a4cb-7e8c46fbdf22
# ╟─6addb17b-1a15-4e66-98d7-03e287664c34
# ╟─593a177f-bf0b-4b05-9408-745793ab2536
# ╟─3f38e8ca-2286-427f-8816-3b8b6cc78c74
# ╟─596e2d59-203c-4e69-985b-f8a82624ef6c
# ╠═465e92c8-1004-47d6-ac4e-69172afad2b0
# ╟─e8d619bc-e37d-437a-96fb-0995aed2f823
# ╟─9dd26db3-e443-46f3-8e18-21eb37b4d5b6
# ╟─36bc055e-3b5b-41a0-863a-53b78a6328d9
# ╟─30c618b6-bd58-4335-bc55-23c16317011d
