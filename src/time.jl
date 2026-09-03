
using Dates: Dates, DateTime, Date
using AstroTime: AstroTime

export mjd
"""
    mjd("2020-01-01")

The modified julian day of a calendar timestamp, for writing epochs down
readably: `mjd("2020-01-01") == 58849.0`.

The timestamp is read on the **TT** (Terrestrial Time) scale, which is the
scale epochs are wanted on — TT and TDB differ by periodic terms below 2 ms,
which no orbit fit can see. The Conventions page of the manual has the full
statement.

!!! warning "This performs no timescale conversion"
    The string is *interpreted* as TT, not converted to it. Neither of the two
    things that separate a raw UTC timestamp at a telescope from a barycentric
    dynamical epoch is applied here:

    * the UTC-to-dynamical offset (leap seconds plus TAI–TT, about 69 s
      today), and
    * the light-travel time from the observer to the solar-system barycentre,
      periodic over the year and reaching ±8.3 minutes.

    So this is the right way to write down a date *you* chose — a plot range,
    a reference epoch, a simulated observation — and the wrong way to convert
    a measured observation timestamp. Measured epochs should reach
    PlanetOrbits already reduced to `BJD_TDB`, which is what instrument
    pipelines deliver.
"""
function mjd(timestamp::AbstractString)
    return timestamp |> 
        AstroTime.TTEpoch |> # Switched from UTC to Terrestrial Time epoch
        AstroTime.modified_julian |>
        AstroTime.days |>
        AstroTime.value;
end
"""
    mjd(Date("2020-01-01"))

The modified julian day of a `Date` or `DateTime`, read on the TT scale. Same
caveats as the string method above: no timescale conversion is performed.
"""
function mjd(date_or_datetime::Union{Date,DateTime})
    return date_or_datetime |> 
        AstroTime.TTEpoch |> # Switched from UTC to Terrestrial Time epoch
        AstroTime.modified_julian |>
        AstroTime.days |>
        AstroTime.value;
end

"""
    mjd()

The modified julian day of right now, on the TT scale.
"""
function mjd()
    return Dates.now() |> 
    AstroTime.TTEpoch |>
    AstroTime.modified_julian |>
    AstroTime.days |>
    AstroTime.value;
end
export mjd


"""
    years2mjd()

Convert from decimal years (e.g. 1995.25) into modified
julian date, rounded to closest second
"""
function years2mjd(decimal_years)
    yr_floor = floor(decimal_years)
    yr_obj = Dates.Date(yr_floor,1,1)
    days = (decimal_years - yr_floor) * Dates.daysinyear(yr_obj)
    days_floor = floor(days)
    ep = AstroTime.TTEpoch(
        Dates.DateTime(yr_floor) + Dates.Day(days_floor) + Dates.Second(round((days-days_floor)*60*60*24))
    )
    return AstroTime.value(AstroTime.modified_julian(ep))
end
export years2mjd

"""
    mjd2date(modified_julian)
    
Get a Date value from a modfied julian day, rounded to closest day

## Examples
```julia
julia> mjd2date(59160.8)
2020-11-08
```
"""
function mjd2date(days)
    return DateTime(
        Dates.DateTime("1858-11-17") + 
        Dates.Day(floor(days)) +
        Dates.Second(round((days-floor(days))*60*60*24))
    )
end
export mjd2date