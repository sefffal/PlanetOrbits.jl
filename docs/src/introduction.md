# Introduction

```@setup 1
# Set a plot theme so that the plots appear nicely on dark or light documenter themes
using Plots
theme(:default;
        framestyle=:box,
    )
```

This package is structured around a representation of an orbit ([`PlanetOrbits.AbstractOrbit`](@ref), and a representation of a "solved" orbit ([`PlanetOrbits.AbstractOrbitSolution`](@ref)).

You start by creating an orbit with known information, e.g. the semi-major axis and eccentricity. You can then query information from this orbit, like its orbital period, mean motion, or periastron (closest approach). Then, you can "solve" the orbit one more times for a given time, eccentric anomaly, or true anomaly.

Let's see how this works.

```@example 1
using PlanetOrbits, Plots
orb = orbit(
    a=1.0, # semi major axis (AU)
    M=1.0, # primary mass (solar masses)
)
```

The [`orbit`](@ref) function accepts many combinations of orbital parameters and returns a subtype of [`PlanetOrbits.AbstractOrbit`](@ref).

We can now query some basic properties about the orbit:
```@example 1
period(orb) # orbital period (days)
```
```@example 1
meanmotion(orb) # Mean motion (radians/yr)
```
```@example 1
periastron(orb) # Epoch of periastron passage (MJD)
```
```@example 1
semiamplitude(orb) # radial velocity semi-amplitude (m/s)
```

We can plot the orbit (more on this in [Plotting](@ref)):
```@example 1
plot(orb)
```

And we can solve the orbit for a given orbital location
```@example 1
sol = orbitsolve_ν(orb, 0.1)        # true anomaly (radians)
sol = orbitsolve_eccanom(orb, 0.1)  # eccentric anomaly (radians)
sol = orbitsolve_meananom(orb, 0.1) # mean anomaly (radians)
```

When constructing an orbit, the location of the planet along its orbit can be specified by `τ`. This is a unitless value between 0 and 1 that represents the fraction of the orbit completed at reference epoch `tref` which by default is the MJD `58849.0`, or 2020-01-01.
```@example 1
orb = orbit(
    a=1.0, # semi major axis (AU)
    M=1.0, # primary mass (solar masses)
    τ=0.2,
);
```

We can now meaningfully solve the location of the planet at a specific time:
```@example 1
t = mjd("2020-07-15") # date as in modified julian days.
sol = orbitsolve(orb, t) # can optionally pass `tref=...` 
```

We can query specifics at this solution:
```@example 1
trueanom(sol) # true anomaly (radians)
```
```@example 1
eccanom(sol) # eccentric anomaly (radians)
```

```@example 1
plot(sol) # defaults to kind=:radvel for RadialVelocityOrbit
```
Notice that we now see a marker at the location found by [`orbitsolve`](@ref).

We can create an orbit with some eccentricity. If not specified, eccentricity and the argument or periapsis default to 0 for any orbit type.
```@example 1
orb = orbit(
    a=1.0, # semi major axis (AU)
    M=1.0, # primary mass (solar masses)
    τ=0.2,
    # New:
    e=0.6, # eccentricity
    ω=2.5, # argument of periapsis (radians)
)
plot(orb) # defaults to kind=:radvel for RadialVelocityOrbit
```

!!! warning "ω convention"
    The convention used in this package is that ω, the argument of periapsis, refers to the **secondary** body. This is in contrast to the typical standard adopted in the radial velocity literature where ω refers to the primary. You can convert by adding or subtracting 180°.

Since we only provided very minimal information to the `orbit` function, we've been receiving a [`RadialVelocityOrbit`](@ref). This object contains sufficient information to calculate the above radial velocity plots, orbital period, etc., but not the 3D position in space.

Let's create a new orbit with a specified inclination and longitude of ascending node.
```@example 1
orb = orbit(
    a=1.0, # semi major axis (AU)
    M=1.0, # primary mass (solar masses)
    τ=0.2,
    e=0.6, # eccentricity
    ω=2.5, # argument of periapsis (radians)
    # New:
    i=0.6, # inclination (radians)
    Ω=2.3, # inclination (radians)
)
```

This time we received a full [`KepOrbit`](@ref). This has the necessary information to solve
the orbit in 2/3D.
```@example 1
plot(orb) # defaults to kind=(:x,:y) for KepOrbit
```
```@example 1
plot(orb, kind=(:x,:y,:z))
```

!!! warning "Cartesian convention"
    The convention used in this package is that x increases to the left (just like right-ascension), and the z increases away from the observer.
    This means all results are ready to be compared to observations where these conventions are typically used; however, it does mean that the space does not follow the standard right-hand rule.
    
We can solve for a time or location as usual.
```@example 1
sol = orbitsolve(orb, mjd("2025-01"))
```
```@example 1
eccanom(sol) # eccentric anomaly (radians)
```

We can also query the cartesian position of the planet in AU:
```@example 1
PlanetOrbits.posx(sol)
```
```@example 1
PlanetOrbits.posy(sol)
```
```@example 1
PlanetOrbits.posy(sol)
```

```@example 1
plot(sol)
```
```@example 1
plot(sol, kind=:x)
```



We can still of course calculate the radial velocity as well.
```@example 1
radvel(sol)
```
```@example 1
plot(sol, kind=:radvel)
```

Finally, we'll specify the parallax distance to the system. This will allow us to plot orbits with angular units as they would appear in the sky from the Earth.
```@example 1
orb = orbit(
    a=1.0, # semi major axis (AU)
    M=1.0, # primary mass (solar masses)
    τ=0.2,
    e=0.6, # eccentricity
    ω=2.5, # argument of periapsis (radians)
    i=0.6, # inclination (radians)
    Ω=2.3, # inclination (radians)
    # New:
    plx=100.0 # parallax distance (milliarcseconds)
)
```

```@example 1
sol = orbitsolve(orb, 0.0)
plot(sol)
```

```
posangle(sol) # position angle offset from barycentre (milliarcseconds)
```

```@example 1
projectedseparation(sol) # separation from barycentre (milliarcseconds)
```

```@example 1
raoff(sol) # right ascension offset from barycentre (milliarcseconds)
```

```@example 1
decoff(sol) # declination offset from barycentre (milliarcseconds)
```

```@example 1
raoff(sol) # right ascension offset from barycentre (milliarcseconds)
```

```@example 1
decoff(sol) # declination offset from barycentre (milliarcseconds)
```

```@example 1
pmra(sol) # instantaneous right ascension velocity from barycentre (milliarcseconds/year)
```

```@example 1
pmdec(sol) # instantaneous declination velocity from barycentre (milliarcseconds/year)
```

```@example 1
accra(sol) # instantaneous right ascension acceleration from barycentre (milliarcseconds/year^2)
```

```@example 1
accdec(sol) # instantaneous declination acceleration from barycentre (milliarcseconds/year^2)
```

## Performance
The [`orbit`](@ref) function is a convenience only for interactive use. It is inefficient since it is not type-stable. Instead, one should use one of the orbit constructors directly.
For example, instead of 
```julia
orb = orbit(
    a=1.0, # semi major axis (AU)
    M=1.0, # primary mass (solar masses)
    τ=0.2,
    e=0.6, # eccentricity
    ω=2.5, # argument of periapsis (radians)
    i=0.6, # inclination (radians)
    Ω=2.3, # inclination (radians)
    plx=100.0 # parallax distance (milliarcseconds)
) # Not type stable
```

Use:
```julia
orb = VisualOrbit(
    a=1.0, # semi major axis (AU)
    M=1.0, # primary mass (solar masses)
    τ=0.2,
    e=0.6, # eccentricity
    ω=2.5, # argument of periapsis (radians)
    i=0.6, # inclination (radians)
    Ω=2.3, # inclination (radians)
    plx=100.0 # parallax distance (milliarcseconds)
) # Type stable
```
This will also suppress the log message.


## Convenience
All functions described above that apply to orbit solutions can be called directly on an orbit along with a time in days:

```@example 1
orb = orbit(
    a=1.0, # semi major axis (AU)
    M=1.0, # primary mass (solar masses)
)
radvel(orb, mjd("2025-01-01"))
```

If you need to calculate many different properties, e.g. both x and y position at a given time/location, it is more efficient to calculate the orbit solution a single time.

## Host calculations
The above calculations treat the planet as a test particle and calculate their displacement/velocity/etc. compared to the two-particle system's barycentre. If you wish to calculate the same properties for the host object, you can additionally supply the mass of the planet.


```@example 1
orb = orbit(
    a=1.0, # semi major axis (AU)
    M=1.0, # primary mass (solar masses)
    i=0.5,
    Ω=2.5,
    plx=100.0
)
sol = orbitsolve(orb, mjd("2025-01-01"))
# Secondary radial velocity
radvel(sol)
```

```@example 1
# Primary radial velocity
radvel(sol, 0.1) # 0.1 solar mass secondary
```

The following show pairs of results for the secondary and the primary:
```@example 1
PlanetOrbits.posx(sol), PlanetOrbits.posx(sol, 0.1)
```

```@example 1
radvel(sol), radvel(sol, 0.1)
```

```@example 1
raoff(sol), raoff(sol, 0.1)
```

```@example 1
accra(sol), accra(sol, 0.1)
```

```@example 1
projectedseparation(sol), projectedseparation(sol, 0.1)
```
```@example 1
posangle(sol), posangle(sol, 0.1)
```
