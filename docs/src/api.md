# API Documentation

The following tables show what functions are supported for what kind of orbit. If you're not sure yet what kind of orbit to use, just use the [`orbit`](@ref) function!
✅ indiciates that a function is available, and ❌ indicates it is not due to the orbit not storing sufficient information. ⚠️ indicates that it could be supoprted, but is not yet implemented.


## Required Parameters
The following table specifies what properties are required to construct each orbit type. Based on this information, different orbit types have different capabilities (described in following tables).

| property  | meaning             | KepOrbit  | Visual{KepOrbit}  | AbsoluteVisual{KepOrbit}  | ThieleInnesOrbit  | RadialVelocityOrbit   | CartesianOrbit    | Visual{CartesianOrbit} |
|---------- | ------------------- |---------- |------------------ |----------------------- |------------------ |---------------------  |----------------   |------------------------|
| M         |                     | ✔️         | ✔️                 | ✔️                      | ✔️                 | ✔️                     | ✔️                 | ✔️                      |
| plx       |                     |           | ✔️                 | ✔️                      | ✔️                 |                       |                   | ✔️                      |
| tp        |                     | ✔️         | ✔️                 | ✔️                      | ✔️                 | ✔️                     | ✔️                 | ✔️                      |
| tref      |                     | ✔️         | ✔️                 | ✔️                      | ✔️                 | ✔️                     | ✔️                 | ✔️                      |
| e         |                     | ✔️         | ✔️                 | ✔️                      | ✔️                 | ✔️                     |                   |                        |
| i         |                     | ✔️         | ✔️                 | ✔️                      |                   |                       |                   |                        |
| ω         |                     | ✔️         | ✔️                 | ✔️                      |                   | ✔️                     |                   |                        |
| Ω         |                     | ✔️         | ✔️                 | ✔️                      |                   |                       |                   |                        |
| A         |                     |           |                   |                        | ✔️                 |                       |                   |                        |
| B         |                     |           |                   |                        | ✔️                 |                       |                   |                        |
| F         |                     |           |                   |                        | ✔️                 |                       |                   |                        |
| G         |                     |           |                   |                        | ✔️                 |                       |                   |                        |
| x         |                     |           |                   |                        |                   |                       | ✔️                 | ✔️                      |
| y         |                     |           |                   |                        |                   |                       | ✔️                 | ✔️                      |
| z         |                     |           |                   |                        |                   |                       | ✔️                 | ✔️                      |
| vx        |                     |           |                   |                        |                   |                       | ✔️                 | ✔️                      |
| vy        |                     |           |                   |                        |                   |                       | ✔️                 | ✔️                      |
| vz        |                     |           |                   |                        |                   |                       | ✔️                 | ✔️                      |
| ref_epoch |                     |           |                   |                        |                   |                       |                   |                        |
| ra        |                     |           |                   |                        |                   |                       |                   |                        |
| dec       |                     |           |                   |                        |                   |                       |                   |                        |
| rv        |                     |           |                   |                        |                   |                       |                   |                        |
| pmra      |                     |           |                   |                        |                   |                       |                   |                        |
| pmdec     |                     |           |                   |                        |                   |                       |                   |                        |

## Properties of Orbits
You can use these functions like `totalmass(orbit)`.

| Function                  | KepOrbit  | Visual{KepOrbit}  | ThieleInnesOrbit  | RadialVelocityOrbit   | CartesianOrbit    | Visual{CartesianOrbit}    |
|----------                 |---------- |------------------ |------------------ |---------------------  |----------------   |------------------------   |
| [`totalmass`](@ref)        | ✅         | ✅                 | ✅                 | ✅                     | ✅                 | ✅                         |
| [`period`](@ref)          | ✅         | ✅                 | ✅                 | ✅                     | ✅                 | ✅                         |
| [`distance`](@ref)        | ❌         | ✅                 | ✅                 | ❌                     | ❌                 | ✅                         |
| [`meanmotion`](@ref)      | ✅         | ✅                 | ✅                 | ✅                     | ✅                 | ✅                         |
| [`eccentricity`](@ref)    | ✅         | ✅                 | ✅                 | ✅                     | ✅                 | ✅                         |
| [`inclination`](@ref)     | ✅         | ✅                 | ✅                 | ❌                     | ✅                 | ✅                         |
| [`semimajoraxis`](@ref)   | ✅         | ✅                 | ✅                 | ✅                     | ✅                 | ✅                         |
| [`periastron`](@ref)      | ✅         | ✅                 | ✅                 | ✅                     | ✅                 | ✅                         |
| [`semiamplitude`](@ref)   | ✅         | ✅                 | ⚠️                 | ✅                     | ✅                 | ✅                         |

## Properties of Orbit Solutions
You can use these functions like `sol = orbitsolve(orbit,mjd("2020-01")); posx(sol)` or  `posx(orbit, mjd("2020-01"))`.

| Function                    | KepOrbit  | Visual{KepOrbit}  | ThieleInnesOrbit  | RadialVelocityOrbit   | CartesianOrbit    | Visual{CartesianOrbit}    |
|----------                   |---------- |------------------ |------------------ |---------------------  |----------------   |------------------------   |
| [`meananom`](@ref)          | ✅         | ✅                 | ✅                 | ✅                     | ✅                 | ✅                         |
| [`trueanom`](@ref)          | ✅         | ✅                 | ✅                 | ✅                     | ✅                 | ✅                         |
| [`eccanom`](@ref)           | ✅         | ✅                 | ✅                 | ✅                     | ✅                 | ✅                         |
| [`posx`](@ref)              | ✅         | ✅                 | ✅                 | ❌                     | ✅                 | ✅                         |
| [`posy`](@ref)              | ✅         | ✅                 | ✅                 | ❌                     | ✅                 | ✅                         |
| [`posz`](@ref)              | ✅         | ✅                 | ✅                 | ❌                     | ✅                 | ✅                         |
| [`velx`](@ref)              | ✅         | ✅                 | ✅                 | ❌                     | ✅                 | ✅                         |
| [`vely`](@ref)              | ✅         | ✅                 | ✅                 | ❌                     | ✅                 | ✅                         |
| [`velz`](@ref)              | ✅         | ✅                 | ✅                 | ✅                     | ✅                 | ✅                         |
| [`raoff`](@ref)             | ❌         | ✅                 | ✅                 | ❌                     | ❌                 | ✅                         |
| [`decoff`](@ref)            | ❌         | ✅                 | ✅                 | ❌                     | ❌                 | ✅                         |
| [`radvel`](@ref)            | ✅         | ✅                 | ✅                 | ✅                     | ✅                 | ✅                         |
| [`posangle`](@ref)          | ✅         | ✅                 | ✅                 | ❌                     | ✅                 | ✅                         |
| [`pmra`](@ref)              | ❌         | ✅                 | ✅                 | ❌                     | ❌                 | ✅                         |
| [`pmdec`](@ref)             | ❌         | ✅                 | ✅                 | ❌                     | ❌                 | ✅                         |
| [`accra`](@ref)             | ❌         | ✅                 | ❌                 | ❌                     | ❌                 | ❌                         |
| [`accdec`](@ref)            | ❌         | ✅                 | ❌                 | ❌                     | ❌                 | ❌                         |
            

## Documentation
```@autodocs
Modules = [PlanetOrbits]
```
