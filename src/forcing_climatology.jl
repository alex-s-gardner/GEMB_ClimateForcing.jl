"""
    forcing_climatology(ds::DimStack)
    forcing_climatology(ds::DimStack, datetime_range::Tuple{DateTime,DateTime})

Compute climatological average forcing from a DimStack.

Creates a single-year average forcing by:
1. Optionally subsetting to `datetime_range`
2. Eliminating leap days (day 366)
3. Excluding partial years
4. Averaging across complete years

Returns a new DimStack with one year of climatological forcing.
Stack metadata is carried forward unchanged so a subsequent
`GEMB.ClimateForcing(ds)` call still finds all required keys.

Matches MATLAB's `forcing_climatology.m`.
"""
function forcing_climatology(ds::DimStack, datetime_range::Tuple{DateTime,DateTime})
    # Subset to the specified date range
    times = collect(lookup(dims(ds, Ti)))
    keep = (times .>= datetime_range[1]) .& (times .<= datetime_range[2])

    tkeep = Ti(times[keep])
    ds_subset = DimStack(
        NamedTuple(k => DimArray(parent(ds[k])[keep], (tkeep,); metadata=metadata(ds[k]))
                   for k in keys(ds));
        metadata = metadata(ds)
    )
    return forcing_climatology(ds_subset)
end

function forcing_climatology(ds::DimStack)
    times = collect(lookup(dims(ds, Ti)))

    # Eliminate leap days (day 366 of the year)
    non_leap = [Dates.dayofyear(t) != 366 for t in times]
    times_noleap = times[non_leap]

    # Count timesteps per year
    years_all = Dates.year.(times_noleap)
    unique_years = sort(unique(years_all))
    counts_per_year = [count(==(yr), years_all) for yr in unique_years]

    # Find years with the maximum number of entries (complete years)
    max_count = maximum(counts_per_year)
    complete_mask = counts_per_year .== max_count
    complete_years = unique_years[complete_mask]

    # Get indices for complete years only
    forcing_index = [yr in complete_years for yr in years_all]
    n_complete_years = length(complete_years)
    steps_per_year = max_count

    times_complete = times_noleap[forcing_index]

    # Reshape into (steps_per_year x n_years) and average
    reshape_avg(arr) = vec(Statistics.mean(reshape(arr[non_leap][forcing_index], steps_per_year, n_complete_years), dims=2))

    # Use times from first complete year
    clim_times = times_complete[1:steps_per_year]
    tdim = Ti(clim_times)

    clim_layers = NamedTuple(
        k => DimArray(reshape_avg(parent(ds[k])), (tdim,); metadata=metadata(ds[k]))
        for k in keys(ds)
    )

    return DimStack(clim_layers; metadata = metadata(ds))
end
