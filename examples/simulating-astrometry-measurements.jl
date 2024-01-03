using PlanetOrbits, Distributions, DataFrames

# The template orbit
template = orbit(
    M = 1.5,
    a = 7.6,
    i = deg2rad(30),
    e = 0.6,
    Ω = 0.7,
    ω = 1.4,
    tp = 0.,
    plx = 100.0,
)

# Epochs specified manually:
epoch = [
    mjd("2020-01-01"),
    mjd("2021-01-01"),
    mjd("2022-01-01"),
]

# Or a range of epochs:
epoch = range(start=mjd("2020-01-01"),step=365,length=3)


astrom = DataFrame(;
    epoch,
    ra=raoff.(template, epoch),
    dec=decoff.(template, epoch),

    # Or:
    # pa=posangle.(template, epoch),# .+ σ_pa .* randn.(),
    # sep=projectedseparation.(template, epoch) .+ σ_sep .* randn.(),
)

# The above will display automatically at a REPL and you can copy the values

## Optional: save to file
using CSV
CSV.write("astrometry.csv", astrom)
