
###################################################
# Quick guide of the package usage by the author
###################################################


#
# Step 1 - Ingest some points
#

# I get given points in coordinate system defined by an SRID (https://en.wikipedia.org/wiki/SRID)
raw_data = [
  573105.43200000003	6086900.3839999996	277.42700000000002	"EPSG"  32755;
  572734.61300000001	6087631.648	        273.80000000000001	"EPSG"  32755;
  572688.821	        6086535.8650000002	267.334	            "EPSG"  32755;
  573667.63899999997	6086370.341			287.68799999999999	"EPSG"  32755
]

# build an srid type from the info in the raw data
auth = symbol(raw_data[1, end-1])
code = raw_data[1, end]
srid = SRID{auth, code}  



# create a vector of the points in this SRID
# Transformation for the SRID_Pos type is handled by Proj4 so anything supported by Proj4 is fine,
# but its up to the user to correctly interpret the x,y,z fields of the SRID_Pos type as they could mean anything (lon lat ... / east north ... / x y ... etc)
Xin = SRID_Pos{srid}(raw_data[:, 1:3], Val{:row})



#
# Step 2 - Transform them into an longitude / latitude / altitude format with the WGS84 ellipsoid (i.e. GPS coordinates) 
#

Xwgs84 = transform(LLA{WGS84}, Xin)




#
# Step 3 - transform them into a local (East, North, Up) coordindate frame
#



# use the center of the points as the center of the local frame
bbox = Bounds(Xwgs84)
lla_ref = center(bbox)


# option a: embed the local frame's origin in the point type.  This saves passing lla_ref between functions
Xlocal_a = transform(ENU{lla_ref}, Xwgs84)
typeof(Xlocal_a)


# option b: don't embed the local frame's origin in the point type.  Use this approach if you use loads of reference points (Julia's type system can only handle so many types)
Xlocal_b = transform(ENU, Xwgs84, lla_ref)
typeof(Xlocal_b)




#
# Step 4 - do something with the points.  For this guide demo we'll just add some new points in between the old ones
#

for pair in combinations(1:length(Xlocal_a), 2)
	push!(Xlocal_a, (Xlocal_a[pair[1]] + Xlocal_a[pair[2]]) / 2.0)
	push!(Xlocal_b, (Xlocal_b[pair[1]] + Xlocal_b[pair[2]]) / 2.0)
end


#
# Step 5 - convert back to a world coordinate frame (ECEF this time to be different)
#

Xecef_a = transform(ECEF, Xlocal_a)           # the reference point is in the type of Xlocal_a
Xecef_b = transform(ECEF, Xlocal_b, lla_ref)  # need to supply the reference point


#
# Step 6 - Lastly go back to the ingest format
#

Xout = transform(SRID_Pos{srid},  Xecef_a)


# N.B. the input srid was actually a utm projection (zone 55) for the WGS84 datum.  We can use this package to select the SRID that best matches another point 
(zone, band) = Geodesy.utm_zone(lla_ref) # I think this is the same for all reference ellipses...
srid_out = Geodesy.utm_srid(lla_ref)    # Geodesy can only do this for the WGS84 datum at present



###############################################
# Warning for the ECEF type
###############################################

# whether ECEF is truly an ECEF point depends on the datum
# e.g.: 
lla_odgb36 = LLA{Geodesy.OSGB36}(0, 0, 0)  # OSGB36 is a good match to the geoid in the UK but not elsewhere
ecef_fake = transform(ECEF, lla_odgb36) 							   # = 6.377563396e6, 0.0, 0.0
# isn't a true ECEF point because the OSGB36 ellipsoid isn't geocentric.

# If in doubt, datum conversions can be done using an SRID_Pos type, 
ecef_srid = transform(SRID(ECEF{WGS84}), lla_odgb36) # = 6.377879171552554e6,-99.12039106890559, 534.423089412207
ecef_wgs84 = convert(ECEF{WGS84}, ecef_srid)         # convert to a type handled by this package (type conversion not a transformation)

# or equivilently
srid = SRID(typeof(lla_odgb36))
lla_srid = convert(SRID_Pos{srid}, lla_odgb36)       # type conversion not a transformation
ecef_wgs84_2 = transform(ECEF{WGS84}, lla_srid)
		
#  Note that the correct SRID won't always be know to this package so the user may need to select the correct one for their datum





