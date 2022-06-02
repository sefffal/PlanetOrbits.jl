
# Symbolic Manipulation
The Symbolics.jl package works fairly well out of the box with PlanetOrbits.jl.
You can create fully or partially symbolic [`KeplerianElements`](@ref) and/or solve for orbits
at a time or true anomaly given by a symbolic `t`.
This could come in use in a few scenarios. For example, if you have an orbit with all parameters known except inclination, you could construct a set of elements with `i` as a symbolic variable.
Solving the orbit using [`orbitsolve`](@ref) would then return a solution with simplified symbolic expressions of `i` that can be evaluated very efficiently for different values.
N.B. this approach is quite messy for a symbolic `e` since Kepler's equation is trancendental.


There is some support for using the Symbolics.jl package. You can create symbolic variables and trace most of the functions defined in this package to get symbolic expressions. 
This is a little slow, and I'm not sure of the applications, but it's neat that it works.

```julia
using Symbolics
@variables t
expr = radvel(elements, t);
```
This works with the KeplerianElements constructors as well if you want to create
a full symbolic set of elements.