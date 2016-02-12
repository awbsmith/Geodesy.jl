### Placing miscellaneuous utm functions here (helpers for common use cases)



# find the utm zone for a specified wgs84 poistion
# the code is stolen from Geograhpiclib::StandardZone()
# http://geographiclib.sourceforge.net/1.19/UTMUPS_8cpp_source.html
# and Geograhpiclib::MGRS::LatitudeBand()
# http://geographiclib.sourceforge.net/2009-02/MGRS_8hpp_source.html
# returns longitude zone and latitude band
function utmzone(lla::Union{LLA, LL})

    # int versions
    ilat = floor(Int64, bound_thetad(lla.lat))
    ilon = floor(Int64, bound_thetad(lla.lon))

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
function utm_srid(lla::Union{LLA{WGS84}, LL{WGS84}})

	# get the band and zone
	(zone, band) = utmzone(lla)
	auth = :EPSG
	code = 32600 + (band > 0 ? 0 : 100) + zone
    return SRID(auth, code)
end


