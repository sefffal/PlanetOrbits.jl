using DirectImages
# using CairoMakie
using GLMakie
using Statistics

using PlanetOrbits

cd(@__DIR__)
img = Float32.(readfits("ABAurigae-SPHERE.fits",1));
# img = centered(img)[-165:165,-165:165]
# img = centered(img)[-250:250,-250:250]
img = centered(img .- mean(img))

imshow2(img, cmap=:magma)
##
° = deg2rad(1)
ot = OrbitalTransformation(
    i = 30°,
    e = 0.0,
    M = 23.0,
    ω = 0.,
    Ω = 60°,
    plx = 2.4642,
    platescale = 1e3/48, # mas/px from scalebar
    dt = 15*365.0
)

using ImageTransformations
i = centered(img)
out = warp(i, ot, axes(i))

imshow2(out, cmap=:magma)
##
heatmap(collect.(axes(out))...,collect(out), axis=(;aspect=DataAspect()),  colormap=:magma,colorrange=(-0.5, 1.2))


##
using Dates
date = Date("2019-12-08")
target_date = Node(Date("2019"))

crop_px_after = 198

img_w = @lift begin

    dt = float(Dates.value($target_date - date))
    ot = OrbitalTransformation(
        i = 30°,
        e = 0.0,
        M = 23.0,
        ω = 0.,
        Ω = 60°,
        plx = 2.4642,
        platescale = 1e3/48, # mas/px from scalebar
        dt = dt
    )
    out = centered(warp(img, ot, axes(img)))
    out[-crop_px_after:crop_px_after,-crop_px_after:crop_px_after]
end

update_theme!(
    font="Arial",
    fontsize=10,
)
fig = Figure(resolution=(720, 720), backgroundcolor=:black)
ax, pl = heatmap(
    fig[1,1],
    img_w,
    axis=(;
        aspect=DataAspect(),
        autolimitaspect = 1,
        tellwidth=true
    ),
    colormap=:magma,
    # :bone_1
    # :seaborn_icefire_gradient
    # colormap=:seaborn_rocket_gradient,
    # :sunset
)


scatter!(
    ax,
    [crop_px_after+1],
    [crop_px_after+1],
    marker="⋆",
    markersize=26,
    color=:white,
    strokewidth=0,
)

colgap!(fig.layout,0)
rowgap!(fig.layout,0)
hidespines!(ax)
hidedecorations!(ax)

# annotations!(ax, "HR8799 from Keck", color=:white, position=(10,3), textsize=40)
# annotations!(ax, "W. Thompson & C. Marois", color=:white, position=(2crop_px_after-135,3), textsize=40,)

annotations!(ax, "AB Aurigae\nSPHERE 1.6 Mm", color=:white, position=(10,3), textsize=32)
annotations!(
    ax,
    @lift(
        Dates.format($target_date, "Y")
    ),
    color=:white,
    position=(10,2crop_px_after-30),
    textsize=32,
)
# annotations!(ax, "W. Thompson\nand C. Marois", color=:white, position=(2crop_px_after-155,3), textsize=32,)



fig

##
for td in Date("2019"):Month(12):Date("2060")
    target_date[] = td
    sleep(0.1)
end
##
record(
    fig,
    "ABAurigae-animation.mp4",
    [
        Date("2019"):Month(12):Date("2100");
        Date("2019"):Month(12):Date("2100");
        Date("2019"):Month(12):Date("2100");
        Date("2019"):Month(12):Date("2100");
    ];
    framerate = 30
) do td
    target_date[] = td
end