
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

Convert from decimal years (e.g. 1995.25) into modified
julian date, rounded to closest second
"""
function years2mjd(decimal_years)
    yr_floor = floor(decimal_years)
    yr_obj = Dates.Year(yr_floor)
    days = (decimal_years - yr_floor) * Dates.daysinyear(yr_obj)
    days_floor = floor(days)
    ep = AstroTime.TTEpoch(
        Dates.DateTime(yr_floor) + Dates.Day(days_floor) + Dates.Second(round((days-days_floor)*60*60*24))
    )
    return AstroTime.value(AstroTime.modified_julian(ep)) # TODO: always days, or do I have to check/convert?
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