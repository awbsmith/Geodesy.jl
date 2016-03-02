
###################################################
# Quick guide for the Geodesy package
###################################################


#
# Step 1 - Grab some points
#

# Starting with points in a coordinate system defined by an SRID (https://en.wikipedia.org/wiki/SRID)
raw_data = [
        573105.43200000003    6086900.3839999996    277.42700000000002    "EPSG"    32755;
        572734.61300000001    6087631.6480000000    273.80000000000001    "EPSG"    32755;
        572688.82100000000    6086535.8650000002    267.334               "EPSG"    32755;
        573667.63899999997    6086370.3410000000    287.68799999999999    "EPSG"    32755
]

# build an srid type from the info in the raw data
auth = symbol(raw_data[1, end-1])
code = raw_data[1, end]
srid = SRID{auth, code}
typealias Utm55S CRS{srid}  # use a type alias if it suits



# create a vector of the points in this SRID
# Transformations for the CRS type are handled by Proj4 package so anything supported by Proj4 is fine, BUT
# its up to the user to correctly interpret the x,y,z fields of the CRS type as they could mean anything (lon lat ... / east north ... / x y ... etc)

Xin = convert(Vector{Utm55S}, raw_data[:, 1:3]; row=true) # = convert(Vector{CRS{srid}}, raw_data[:, 1:3]; row=true)



#
# Step 2 - Transform the points into an longitude / latitude / altitude format with the WGS84 ellipsoid (i.e. GPS coordinates)
#
# N.B.  The below step requires Geodesy to recognise the Proj4 projection for type LLA{WGS84}.  If you're using a point type
# that geodesy doesn't know the Proj4 projection for, overload Geodesy.get_projection() to return the appropriate Proj4 projection e.g.:
#    Geodesy.get_projection(::Type{LLA{WGS84}}) = Proj4.Projection("+proj=longlat +datum=WGS84 +no_defs")
#

Xwgs84 = geotransform(LLA{WGS84}, Xin)




#
# Step 3 - transform them into a local (East, North, Up) coordindate frame
#

# use the center of the points as the center of the local frame
bbox = Bounds(Xwgs84)
lla_ref = center(bbox)


# option a: embed the local frame's origin in the point type.  This saves passing around lla_ref
Xlocal_a = geotransform(ENU{lla_ref}, Xwgs84)
println(typeof(Xlocal_a))


# option b: don't embed the local frame's origin in the point type.  Use this approach if you use loads of reference points 
Xlocal_b = geotransform(ENU, Xwgs84, lla_ref)
println(typeof(Xlocal_b))




#
# Step 4 - do something with the points.  For this guide we'll just add some new points in between the old ones
#

for pair in combinations(1:length(Xlocal_a), 2)
    push!(Xlocal_a, (Xlocal_a[pair[1]] + Xlocal_a[pair[2]]) / 2.0)
    push!(Xlocal_b, (Xlocal_b[pair[1]] + Xlocal_b[pair[2]]) / 2.0)
end


#
# Step 5 - convert back to a world coordinate frame (ECEF this time to be different)
#

Xecef_a = geotransform(ECEF, Xlocal_a)           # the reference point is in the type of Xlocal_a
Xecef_b = geotransform(ECEF, Xlocal_b, lla_ref)  # need to supply the reference point


#
# Step 6 - Lastly go back to the original format
#

Xout = geotransform(CRS{srid},  Xecef_a)


# N.B. the input srid was actually a utm projection (zone 55) for the WGS84 datum.  We can use this package to select the SRID that best matches another point 
(zone, band) = Geodesy.utm_zone(lla_ref) # I think this is the same for all reference ellipses...
srid_out = Geodesy.utm_srid(lla_ref)     # Geodesy can only do this for the WGS84 datum at present



###############################################
# ECEF Type Warning 
###############################################

# whether ECEF is truly an ECEF point depends on the datum
# e.g.: 
lla_osgb36 = LLA{Geodesy.OSGB36}(0, 0, 0)          # OSGB36 is a good match to the Earth in the UK but not elsewhere
ecef_fake = geotransform(ECEF, lla_osgb36)         # = 6.377563396e6, 0.0, 0.0
# isn't a true ECEF point because the OSGB36 ellipsoid isn't geocentric.


# If in doubt, datum conversions can be using the CRS type to get Proj4 to do it, e.g. 
ecef_srid = geotransform(SRID(ECEF{WGS84}), lla_osgb36)     # = 6.377879171552554e6,-99.12039106890559, 534.423089412207
ecef_wgs84 = convert(ECEF{WGS84}, ecef_srid)             # convert to a type native to this package (type conversion not a transformation)

# or equivilently
srid = SRID(lla_odgb36)
lla_srid = convert(CRS{srid}, lla_osgb36)           # type conversion not a transformation
ecef_wgs84_2 = geotransform(ECEF{WGS84}, lla_srid)            # = 6.377879171552554e6,-99.12039106890559, 534.423089412207
        

# geotransform(WGS84, 



