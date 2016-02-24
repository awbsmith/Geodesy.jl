# Geodesy

#[![Build Status](https://travis-ci.org/JuliaGeo/Geodesy.jl.svg?branch=master)](https://travis-ci.org/JuliaGeo/Geodesy.jl)
#[![Coverage Status](http://img.shields.io/coveralls/JuliaGeo/Geodesy.jl.svg)](https://coveralls.io/r/JuliaGeo/Geodesy.jl)

Geodesy is designed to work with world locations various coordinate systems. The code has been split out from [OpenStreetMap.jl](https://github.com/tedsteiner/OpenStreetMap.jl), and functionality expanded.

![Coordinate Reference systems](http://www.crs-geo.eu/SharedDocs/Bilder/CRS/schema-crs-datum-cs,property=default.gif)
[www.crs-geo.eu](http://www.crs-geo.eu/nn_124224/crseu/EN/CRS__Overview/definition-crs__node.html)

The above image gives a quick picture of the components of the coordinate reference systems used in geodesy.  This Geodesy package is intended for use with the "Coordinate System" subtypes, however the [Proj4 package](https://github.com/FugroRoames/Proj4.jl) is used as a backend top allow tranforming to / from full "Coordinate _Reference_ Sytems". 

### "Coordinate Reference System" Types

The below types are parameterised by a spatial reference ID (SRID). Note that's SRID's are used to describe more than just coordinate reference systems, so take care when selecting selected a suitable SRID.

1. `CRS` - The coordinate _reference_ system point type.  This type should be used for operations that require knowledge of the datum, e.g. swapping between datums. Transformations involving this type are perfromed by Proj4 as a full understanding of coordinate reference systems is beyond the scope of this package. It's left to the user to correctly interpret the fields of this type as different SRID use different coordinate systems ([lat, long, height] / [x, y, z] / [false east, false north, height] / etc).


### "Coordinate System" Types

The below types are parameterised by a reference ellipse. Tags for common ellipse's are provided, and custom ellipse's can also be used. Tag types for common coordinate reference systems are also provided as a parameter, but they are only a convenience for getting the coordinate reference system's ellipsoid (e.g. OSGB36 uses the [AIRY ellipsoid](https://en.wikipedia.org/wiki/Ordnance_Survey_National_Grid#General)).

1. `LLA`   - Longitude latitude height coordinate

2. `LL`    - Longitude latitude point coordinate to be on the ellipse's surface

3. `ECEF`  - X Y Z coordinate.  Note that ECEF has a [meaning beyong a Cartesian cooridnate system](https://en.wikipedia.org/wiki/ECEF), however is used in this package to mean only a Cartesian cooridnate system with origin equal to the ellipse's center.  More detail on this is provided below.


The below type is parameterized by an `LL` point which determines the origin and axis directions of the local coordinate system.  The `LL` parameter can be omitted from the template and passed to function seperately if desired.

4. `ENU`   - East North Up poistion in the local coordinate frame.


### Transformations vs Conversions

1. `convert` is used for value preserving transformations, e.g.:
```julia
	convert(CRS{SRID(LLA{WGS84})}, LLA{WGS84}(0.0,0.0,0.0)) 
```
is appropriate as the point value remains the same.  `convert` should not change values!

1. `transform` is used for value modifying transformations, e.g.:
```julia
transform(ECEF, LLA{WGS84}(0.0,0.0,0.0)) # = ECEF{WGS84}(6.378137e6,0.0,0.0)
```


### Examples

##### Convert bewteen "Coordinate System" types

```julia

# start with an LLA point
lla_wgs84 = LLA{WGS84}(.1167, 51.5, 0.0)

# convert to a Cartesian coordinate frame
ecef_wgs84 = transform(ECEF, lla_wgs84)   # or equivilently:
ECEF(lla_wgs84)

# convert to a local coordinate frame
lla_ref = LLA{WGS84}(.1168, 51.5, 0.0)
enu_ref = transform(ENU{lla_ref}, lla_wgs84)    # = ENU{LL{WGS84}(0.1167,51.5)}(6.944051663969915,4.742430141962267e-6,-3.772299267203877e-6)
enu_noref = transform(ENU, lla_ref, lla_wgs84)) # = ENU{???}(6.944051663969915,4.742430141962267e-6,-3.772299267203877e-6)

# convert from the local coordinate frame to the LLA coordinate frame
lla = transform(LLA, enu_ref)
lla = transform(LLA, enu_ref, lla_wgs84)
```



##### Convert a point in a specified SRID to a "Coordinate System" Type

```julia

# create a point in a coordinate reference system
srid = SRID{:EPSG, 32755}    													# WGS84 datum [UTM](https://en.wikipedia.org/wiki/Universal_Transverse_Mercator_coordinate_system) zone 55 South

utm = CRS{srid}(573105.43200000003, 086900.3839999996, 277.42700000000002) 		# encodes a 2D projection of the Earth's surface

# convert to a coordinate sytem type
lla_wgs84 = transform(LLA{WGS84}, utm)  # the SRID corresponding to LLA{WGS84} is known to Geodesy (see known_srids.jl).  Otherwise, 

lla_wgs84 = convert(LLA{WGS84}, transform(SRID{:EPSG, 4326}, utm))  # EPSG4326 is SRID for for the WGS84 LLA coordinate reference system

```


##### A note on the ECEF type

An ECEF coordinate system's origin should be the Earth's center of mass and have axes aligned with the International Reference Pole (IRP) and International Reference Meridian.  The ECEF point type used
in this package is simply a Cartesian coordinate system with an origin set to the ellipse parameters center.  

As an example, the OSGB36 datum is not centered on the Earth's center of mass, and so an ECEF{OSGB36} point is not a true ECEF point.  This can be demonstrated by compare the transformations of coordinate system and coordinate _reference_ system types:

```julia

# OSGB36 is a good match to the Earth in the UK but not elsewhere
import Geodesy: OSGB36
lla_osgb36 = LLA{OSGB36}(0, 0, 0)  		

# converting to ECEF using only the coordinate system info
ecef_fake = transform(ECEF, lla_osgb36) 		# = 6.377563396e6, 0.0, 0.0

# converting to ECEF using the full datum info
srid = SRID(LLA{OSGB36})
crs_osgb36 = convert(CRS{srid}, lla_osgb36)

# There's different estimates of the Earth's center of mass etc, use the WGS84 (GPS) estimate here
ecef_srid = transform(SRID(ECEF{WGS84}), crs_osgb36) 	# = 6.377879171552554e6,-99.12039106890559, 534.423089412207
ecef_wgs84 = convert(ECEF{WGS84}, ecef_srid)         	

# N.B. you can do this too
ecef_wgs84 = convert(ECEF{WGS84}, transform(SRID(ECEF{WGS84}), lla_osgb36))



```



