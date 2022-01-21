
using Dates, AstroTime

export mjd
"""
    mjd("2020-01-01")

Get the modfied julian day of a date, or in general a UTC
timestamp.
"""
function mjd(timestamp::AbstractString)
    return timestamp |> 
        TTEpoch |> # Switched from UTC to Terrestrial Time epoch
        modified_julian |>
        days |>
        value;
end

"""
    mjd()

Get the current modified julian day of right now.
"""
function mjd()
    return Dates.now() |> 
    UTCEpoch |>
    modified_julian |>
    days |>
    value;
end
export mjd


"""
    years2mjd()

Convert from fractional years (e.g. 1995.25) into modified
julian date.
"""
years2mjd(years) = 365.25*(years - 1858 - 321/365.25)
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
    return Dates.Date("1858-11-17") + Dates.Day(round(days))
end
export mjd2date