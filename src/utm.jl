### Placing miscellaneuous utm functions here (helpers for common use cases)



# find the utm zone for a specified wgs84 poistion
# the code is stolen from Geograhpiclib::StandardZone()
# http://geographiclib.sourceforge.net/1.19/UTMUPS_8cpp_source.html
# and Geograhpiclib::MGRS::LatitudeBand()
# http://geographiclib.sourceforge.net/2009-02/MGRS_8hpp_source.html
# returns longitude zone and latitude band
# TODO: is this correct for ellipsoids other than the wgs84 ellipsoid?
function utmzone(LLA_wgs84::Union{LLA{WGS84}, LL{WGS84}})

    # int versions
    ilat = floor(Int64, bound_thetad(LLA_wgs84.lat))
    ilon = floor(Int64, bound_thetad(LLA_wgs84.lon))

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

macro add_type(str)
	return 
end

# find the EPSG SRID for a utm projection of a specified wgs84 poistion
function utm_srid(LLA_wgs84::Union{LLA{WGS84}, LL{WGS84}})

	# get the band and zone
	(zone, band) = utmzone(LLA_wgs84)
    sym = symbol(@sprintf("EPSG32%i%02i", (band > 0 ? 6 : 7), zone))

    return sym
end


