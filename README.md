# Geodesy

#[![Build Status](https://travis-ci.org/JuliaGeo/Geodesy.jl.svg?branch=master)](https://travis-ci.org/JuliaGeo/Geodesy.jl)
#[![Coverage Status](http://img.shields.io/coveralls/JuliaGeo/Geodesy.jl.svg)](https://coveralls.io/r/JuliaGeo/Geodesy.jl)

Geodesy is designed to work with world locations various coordinate systems. 

![Coordinate Reference systems](http://www.crs-geo.eu/SharedDocs/Bilder/CRS/schema-crs-datum-cs,property=default.gif)
[www.crs-geo.eu](http://www.crs-geo.eu/nn_124224/crseu/EN/CRS__Overview/definition-crs__node.html)

The above image gives a quick picture of the components of the coordinate reference systems used in geodesy.  This Geodesy package is intended for use with the "Coordinate System" subtypes, however the [Proj4 package](https://github.com/FugroRoames/Proj4.jl) is used as a backend to allow transforming to / from / between "Coordinate _Reference_ Systems". 

Roadmap note: Its intended to define the coordinate reference system types in backend packages (Proj4 etc), with the Geodesy package provided common interface methods for them (functions to overload etc).  This change would mean Proj4 depends on Geodesy instead of Geodesy depending on Proj4.

### "Coordinate Reference System" Types

The below type is parameterised by a [spatial reference ID (SRID)](https://en.wikipedia.org/wiki/SRID). Note that SRIDs are used to describe more than just coordinate reference systems, so take care when selecting selected a suitable SRID.

1. `CRS` - The coordinate _reference_ system point type.  This type should be used for operations that require knowledge of the datum, e.g. swapping between datums. Transformations involving this type are perfromed by Proj4 as a full understanding of coordinate reference systems is beyond the scope of this package. It's left to the user to correctly interpret the fields of this type as different SRIDs use different coordinate systems ([lat, long, height] / [x, y, z] / [false east, false north, height] / etc).


### "Coordinate System" Types

The below types are parameterised by a reference datum. Note that the coordinate system types only understands the datum's ellipse; however using a datum type as a parameter is a convenient way to get the reference ellipse while also stopping direct comparison of points from different datums that may use the same ellipse.  A discussion on datums vs ellipses is given later in this readme.

Some common ellipses and datums are provided, and custom ellipse's can also be used. For a list of all predefined datums, use `Geodesy.get_datums()`

1. `LLA`   - Longitude latitude height coordinate

2. `LL`    - Longitude latitude point coordinate to be on the ellipse's surface

3. `ECEF`  - X Y Z coordinate.  Note that ECEF has a [meaning beyong a Cartesian coordinate system](https://en.wikipedia.org/wiki/ECEF), however it is used in this package to mean only a Cartesian coordinate system with origin equal to the ellipse's center.  This is a consequence of this package only using the ellipse part of the datum.  An example of this is given below.


The below type is parameterized by an `LL` point which determines the origin and axis directions of the local coordinate system.  The `LL` parameter can be omitted from the template and passed to function seperately if desired.

4. `ENU`   - East North Up position in the local coordinate frame.


### Transformations vs Conversions

1. `convert` is used for value preserving transformations, e.g. change the position's type to include the coordinate reference system:
```julia
	lla_wgs84 = LLA{WGS84}(0.0,0.0,0.0)
	lla_wgs84_srid = SRID(LLA{WGS84})
	convert(CRS{lla_wgs84_srid}, lla_wgs84) # = CRS{EPSG4326}(0.0,0.0,0.0)
```
	convert` will not change values!

2. `geotransform` is used for value modifying transformations, e.g.:
```julia
	lla_wgs84 = LLA{WGS84}(0.0,0.0,0.0)
	geotransform(ECEF, lla_wgs84) # = ECEF{WGS84}(6.378137e6,0.0,0.0)
```


### Examples

##### Convert between "Coordinate System" types

```julia

# start with an LLA point
lla_wgs84 = LLA{WGS84}(.1167, 51.5, 0.0)

# convert to a Cartesian coordinate frame
ecef_wgs84 = geotransform(ECEF, lla_wgs84)   # or equivilently:
ECEF(lla_wgs84)

# convert to a local coordinate frame
lla_ref = LLA{WGS84}(.1168, 51.5, 0.0)  # center the local frame on this point
enu_wref = geotransform(ENU{lla_ref}, lla_wgs84)    # = ENU{LL{WGS84}(0.1168,51.5)}(-6.944051663969915,4.742430141962267e-6,-3.772299267203877e-6)
enu_noref = geotransform(ENU, lla_wgs84, lla_ref)  # = ENU{???}(-6.944051663969915,4.742430141962267e-6,-3.772299267203877e-6)

# convert from the local coordinate frame to the LLA coordinate frame
lla = geotransform(LLA, enu_wref)              # =LLA{WGS84}(.1167, 51.5, 0.0)
lla = geotransform(LLA, enu_noref, lla_ref)    # =LLA{WGS84}(.1167, 51.5, 0.0)
```



##### Convert a point in a specified SRID to a "Coordinate System" Type

```julia

# create a point in a coordinate reference system
srid = SRID{:EPSG, 32755}    													# WGS84 datum [UTM](https://en.wikipedia.org/wiki/Universal_Transverse_Mercator_coordinate_system) zone 55 South
utm = CRS{srid}(573105.43200000003, 086900.3839999996, 277.42700000000002) 		# a point in the above zone

# convert to a coordinate sytem type
lla_wgs84 = geotransform(LLA{WGS84}, utm)  							   # the SRID corresponding to LLA{WGS84} is known to Geodesy (see known_srids.jl).  Otherwise, 
lla_wgs84 = convert(LLA{WGS84}, geotransform(SRID{:EPSG, 4326}, utm))  # EPSG4326 is the SRID for for the WGS84 LLA coordinate reference system

```

##### Perform transformations on custom types

```julia

	import Geodesy: geodify

	# define a custom type
	immutable CustomLLA
		time::DateTime
		txt::ASCIIString
		lon::Float64
		lat::Float64
		alt::Float64
	end

	# define its conversion to a geodesy type
	geodify(X::CustomLLA) = LLA{WGS84}(X.lon, X.lat, X.alt)

	# instantiate it
	custom_lla = CustomLLA(now(), "test input", .1167, 51.5, 0.0)
	geodify(custom_lla) # = Geodesy.LLA{WGS84}(0.1167,51.5,0.0)

	# and transform
	ecef = geotransform(ECEF, custom_lla) # = ECEF{WGS84}(3.9786402778542214e6,8103.702688750244,4.968362457291028e6)

	# make the output a custom type as well
	immutable CustomECEF
		time::DateTime
		txt::ASCIIString
		x::Float64
		y::Float64
		z::Float64
	end

	# define the gedoesy type required to construct it
	geodify(::Type{CustomECEF}) = ECEF{WGS84}

	# define a constructor that takes a geodesy type input, as well as the input to the transform 
	import Base.call  
	Base.call(::Type{CustomECEF}, ecef::ECEF, X::CustomLLA) = CustomECEF(X.time, X.txt, ecef[1], ecef[2], ecef[3])

	# and transform
	custom_ecef = geotransform(CustomECEF, custom_lla) # = CustomECEF(2016-02-25T16:13:59,"test input",3.9786402778542214e6,8103.702688750244,4.968362457291028e6)

```


#### A note on datums vs ellipses

Datums contain more information than just the reference ellipse.  Two different datums can use the same ellipse but have differenent centers (relative to a point on the Earth's surface, or have different reference poles / meridians.  Comparing points in the different datums also requires having a transforming to align them.  As an example,


```julia

# The [vn2000 datum](https://epsg.io/4756-5194) for Vietnam, which uses the WGS84 ellipsoid 
vn2000 = SRID{:EPSG, 4756}  # SRID to parameterize the Coordinate Reference System

# define this datum for use with Geodesy
immutable VN2000 <: Geodesy.Datum; end
import Geodesy: ellipsoid
ellipsoid(::Type{VN2000}) = ellipsoid(Geodesy.WGS84_ELLIPSE)  # overload the function to retrieve this datum's ellipsoid

# define an LLA point in this datum
lla_vn2000 = LLA{VN2000}(0, 0, 0)
crs_vn2000 = convert(CRS{vn2000}, lla_vn2000)

# convert this point to the WGS84 datum
lla_wgs84 = geotransform(LLA{WGS84}, crs_vn2000)   # = LLA{WGS84}(-0.0003528546332018209,-0.0010055677147903898,-192.75199438724667)

```

We can check the Proj4 well known text for the VN2000 datum and see the *+towgs84* transformation parameters

```Julia
Geodesy.proj4_str(vn2000) # = "+proj=longlat +ellps=WGS84 +towgs84=-192.873,-39.382,-111.202,-0.00205,-0.0005,0.00335,0.0188 +no_defs"
```


#### A note on the ECEF type

An ECEF coordinate system's origin should be the Earth's center of mass and have axes aligned with the International Reference Pole and International Reference Meridian ([ECEF](https://en.wikipedia.org/wiki/ECEF)].  Since the coordinate system types only use the ellipse part of the datum, they have no information to align the datums axis to the International Reference Meridian etc. 

As a result the ECEF point type used in this package is simply a Cartesian coordinate system with an origin set to the parameterising ellipse's center, with axes specified by the ellipse which could point anywhere on the Earth's surface.  The tag ECEF is choosen for pragmatism.

As an example:

```julia
lla_osgb36 = LLA{Geodesy.OSGB36}(0, 0, 0)  		# OSGB36 is a good match to the Earth in the UK but not elsewhere, and is not centered on the Earth's center of mass (i.e. not a [geocentric datum](http://support.esri.com/en/knowledgebase/GISDictionary/term/geocentric%20datum))
cart = geotransform(ECEF, lla_osgb36) 		            # = 6.377563396e6, 0.0, 0.0

# we can compare the above to the ECEF position using the WGS84 datum (a geocentric datum)
ecef_crs = geotransform(SRID(ECEF{WGS84}), lla_osgb36) 	# = 6.377879171552554e6,-99.12039106890559, 534.423089412207
ecef = convert(ECEF{WGS84}, ecef_crs)         	        # convert to a type native to this package (type conversion not a transformation)

# or equivilently
srid = SRID(lla_osgb36)
lla_crs = convert(CRS{srid}, lla_osgb36)       	        # type conversion not a transformation
geotransform(ECEF{WGS84}, lla_crs)			            # = 6.377879171552554e6,-99.12039106890559, 534.423089412207

# the difference
dist = norm(Vector(cart) - Vector(ecef))                # = 628.6072621385788

# which is why we cant do this directly
cart - ecef # = error


```














