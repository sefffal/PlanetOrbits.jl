
using Dates: Dates, DateTime, Date
using AstroTime: AstroTime

export mjd
"""
    mjd("2020-01-01")

Get the modfied julian day of a date, or in general a UTC
timestamp.
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

Get the modfied julian day of a Date or DateTime object.
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

Get the current modified julian day of right now.
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
    return DateTime(
        Dates.DateTime("1858-11-17") + 
        Dates.Day(round(days)) +
        Dates.Second(round((days-floor(days))*60*60*24))
    )
end
export mjd2date