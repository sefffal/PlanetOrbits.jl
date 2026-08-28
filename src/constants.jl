# ---------------------------------------------------
# Constants
# ---------------------------------------------------

# radians <-> arcseconds
const rad2as = 360*3600/2π
const as2rad = 2π/(360*3600)

# parsecs <-> astronomical units
const pc2au = 360*3600/2π # Exact (but irrational) by IAU definition
const au2pc = 1/pc2au

# astronomical units <-> metres (IAU definition)
const au2m = 1.495978707e11
const m2au = 1/au2m

# days <-> seconds
const day2sec = 86400
const sec2day = 1/day2sec

"""
    year2day_julian

Days per Julian year, 365.25 exactly (IAU). Exported because the `P` orbital
element is in **days** while periods are usually quoted, prior'd and plotted
in years, so user code converts constantly:

    P = P_yr * year2day_julian     # sampled a period in years

Deliberately *not* named `year2day`: this package carries two years that differ
by 1.9e-5 — the Julian year here and
[`kepler_year_to_julian_day_conversion_factor`](@ref) (365.2568983840419 days),
the period of a 1 M⊙ / 1 AU orbit under the IAU nominal constants. Mixing them
produces results that look right and quietly drift. This is the one you want
for a period expressed in years; the other only appears inside Kepler's third
law.
"""
const year2day_julian = 365.2500  # IAU definition of a Julian year

"""
    day2year_julian

Julian years per day, `1/year2day_julian`. See [`year2day_julian`](@ref).
"""
const day2year_julian = 1/year2day_julian

export year2day_julian, day2year_julian

# years <-> seconds. IAU definition, ie using Julian years
const year2sec_julian = year2day_julian*day2sec
const sec2year_julian = 1/year2sec_julian

"""
    PlanetOrbits.kepler_year_to_julian_day_conversion_factor

The "Kepler year" in julian days, 365.2568983840419 — the period of a
hypothetical planet orbiting a 1 M⊙ star at exactly 1 AU under the IAU
nominal constants. This is what appears inside Kepler's third law:

    period_days = √(a³ / M) * kepler_year_to_julian_day_conversion_factor

It exists because the IAU definitions of the AU and the solar mass do not
combine to give an orbital period of exactly one julian year. It therefore
differs from [`year2day_julian`](@ref) (365.25) by 1.9e-5 — small enough that
substituting one for the other looks right and quietly drifts. Use this one
only in Kepler's third law; use `year2day_julian` for a period quoted in
years.

Unexported, and deliberately verbose, so that the distinction has to be made
on purpose.

Derivation, from Gilles Otten (thank you for tracking this down): with
`G*M_sun = 1.3271244e20 m³/s²` and `1 AU = 149_597_870_700 m` by definition,
`sqrt(4π²/1.3271244e20 * 149597870700³) / 86400 = 365.2568983840419` julian
days.
"""
const kepler_year_to_julian_day_conversion_factor = 365.2568983840419 # julian days

# GM for M = 1 M⊙ in AU³ per **julian** year², i.e. 4π² (which is GM per
# *kepler* year²) rescaled to the time unit `Row`'s mean motion `n` and
# velocity factor `J` actually use. The two years differ by only 1.9e-5, so
# mixing them produces results that look right and quietly drift — always
# reach for this rather than a bare 4π².
const GM_sun_au3_julianyr2 =
    4π^2 * (year2day_julian / kepler_year_to_julian_day_conversion_factor)^2

"""
    msun

One IAU solar mass, in the mass unit used throughout PlanetOrbits (which *is*
the solar mass, so `msun == 1.0`). Provided so masses read naturally:
`Body(mass=1.2msun)`.
"""
const msun = 1.0

"""
    mjup

One jupiter mass in solar masses (exact ratio of the IAU nominal values of
GM_jup and GM_sun): `Body(mass=5.3mjup)`.
"""
const mjup = 1.2668653e17/1.3271244e20 # == 0.0009545942339693249

"""
    mearth

One earth mass in solar masses (exact ratio of the IAU nominal values of
GM_earth and GM_sun): `Body(mass=23mearth)`.
"""
const mearth = 3.986004e14/1.3271244e20

export msun, mjup, mearth

# Constants used by the absolute-frame 3D motion compensation.
#
# NB: v0.11 (orbit-absolute.jl) carried `c = 2.998e8` m/s and a separate hardcoded
# `2.99792e5` km/s in the light-travel term — 2.5e-5 and 1.5e-6 relative errors
# respectively. Their net effect was a 1.3e-6 relative error on the barycentric
# light-travel time (0.6 s on Barnard's-star-like 4.6e5 s), almost all of which
# is degenerate with the period. v1 uses the exact IAU value throughout, so
# absolute-frame epochs differ from v0.11 at that level by design.
const c_light_ms = 2.99792458e8  # m/s, exact by SI definition
const pc2km = 3.08567758149137e13
const pc2m = pc2km * 1000
const rad2as_206265 = 180 / π * 60 * 60  # rad -> arcsec factor as written in v0.11
const one_over_pc2km_sec2yr = 1.022712165045694034700736065713114217745793404987068055763763987835564887975633e-06

# Light travel time across one parsec [s]. The single definition every
# light-travel term routes through, so there is exactly one speed of light in
# the package. (v0.11 had two: `2.998e8` m/s and a separately hardcoded
# `2.99792e5` km/s in the light-travel term.)
const pc2sec_light = pc2m / c_light_ms

# Speed of light in the units the propagation frame works in. Per-body
# light-travel-time retardation is a delay of z/c with z in AU, so this is the
# constant that matters for it — it needs no distance and no frame at all.
const c_au_per_julianyr = c_light_ms * year2sec_julian / au2m  # ≈ 63241.08 AU/yr
const c_au_per_day = c_light_ms * day2sec / au2m

# radians -> milliarcseconds, i.e. the numerator of the AU -> mas conversion at
# a given distance in AU: cart2angle = rad2mas / d_au.
const rad2mas = rad2as * 1000
