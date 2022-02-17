
# Units & Conventions

The main constructor, [`KeplerianElements`](@ref), accepts the following parameters:
- `a`: Semi-major axis in astronomical units (AU)
- `i`: Inclination in radians
- `e`: Eccentricity in the range [0, 1)
- `τ`: Epoch of periastron passage, in fraction of orbit [0,1]
- `M`: Graviataion parameter of the central body, expressed in units of Solar mass.
- `ω`: Argument of periastron
- `Ω`: Longitude of the ascending node, radians.
- `plx`: Distance to the system expressed in milliarcseconds of parallax.

Thee parameter `τ` represents the epoch of periastron passage as a  fraction of the planet's orbit between 0 and 1. This follows the same convention as Orbitize! and you can read more about their choice in ther FAQ.

Parameters can either be specified by position or as keyword arguments (but not a mix).

See this PDF for a detailed derivation of projected position, velocity, and acceleration from these coordinates: [Derivation.pdf](assets/orbit_coordinate_notes.pdf)

There is also a convenience constructor [`KeplerianElementsDeg`](@ref) that accepts `i`, `ω`, and `Ω` in units of degrees instead of radians.


```@raw html
<img src="https://docs.exoplanet.codes/en/latest/_images/orbit3D.png" style="width:300px"/>
```
**Orbit Convenctions Schematic. Credit: [exoplanet.py](https://docs.exoplanet.codes/en/latest/).**

This diagram from exoplanet.py is a good reference for the conventions used by this package with one exception: we flip the z-coordinate such that radial velocity is positive increasing away from the Earth.
This does mean that the z-coordinate does not follow the right-hand rule as one might expect.



