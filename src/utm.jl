### Placing miscellaneuous utm functions here (helpers for common use cases)



# find the utm zone for a specified LLA position
# TODO: Make sure this is the same for all Ellipsoids
"""
Find the utm zone and band for the input type.  get_lat() and get_lon() must be defined for the input
"""
function utm_zone(lla)

    # int versions
    ilat = floor(Int64, bound_thetad(get_lat(lla)))
    ilon = floor(Int64, bound_thetad(get_lon(lla)))

    # get the latitude band
    band = max(-10, min(9,  fld((ilat + 80), 8) - 10))

    # and check for weird ones
    zone = fld((ilon + 186), 6);
    if ((band == 7) && (zone == 31) && (ilon >= 3))
        zone = 32
    elseif ((band == 9) && (ilon >= 0) && (ilon < 42))
        zone = 2 * fld((ilon + 183), 12) + 1
    end

    return (zone, band)

end

# TODO: overload with other ellipses as needed
"""
Find the SRID for utm zone that the lat lon point input lies in.  This is only valid for the WGS84 datum.

Returns the authority and the code
"""
function utm_srid(lla::Union{LLA{WGS84}, LL{WGS84}})

    # get the band and zone
    (zone, band) = utm_zone(lla)
    auth = :EPSG
    code = 32600 + (band > 0 ? 0 : 100) + zone
    return (auth, code)
end


