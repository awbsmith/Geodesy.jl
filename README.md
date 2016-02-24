# Geodesy

#[![Build Status](https://travis-ci.org/JuliaGeo/Geodesy.jl.svg?branch=master)](https://travis-ci.org/JuliaGeo/Geodesy.jl)
#[![Coverage Status](http://img.shields.io/coveralls/JuliaGeo/Geodesy.jl.svg)](https://coveralls.io/r/JuliaGeo/Geodesy.jl)

Geodesy is designed to work with world locations various coordinate systems. The code has been split out from [OpenStreetMap.jl](https://github.com/tedsteiner/OpenStreetMap.jl), and functionality expanded.

![Coordinate Reference systems](http://www.crs-geo.eu/SharedDocs/Bilder/CRS/schema-crs-datum-cs,property=default.gif)
[www.crs-geo.eu](http://www.crs-geo.eu/nn_124224/crseu/EN/CRS__Overview/definition-crs__node.html)

The above image gives a quick picture of the components of the coordinate reference systems used in geodesy applications.  This Geodesy package is intended for use with the "Coordinate System" subtypes, however the [Proj4 package](https://github.com/FugroRoames/Proj4.jl) is used as a backend top allow tranforming to / from full "Coordinate _Reference_ Sytems" 

### "Coordinate Reference System" Types

1. `CRS` - The coordinate _reference_ system type, parameterized by a spatial reference ID (SRID).  This type should be used for operations such as finding an equivalent point with a different datum.  Note   that's SRID's are used to desribe more than just coordinate reference systems, so take care when selecting selected a suitable SRID.  Transformations involving this type are perfromed by Proj4 as a full understanding of coordinate reference systems is beyoing the scope of this package.

		It's left to the user to correctly interpret the fields of this type as different SRID use different coordinate systems (lat, long, height / x, y, z / false east, false north, height / etc)


### "Coordinate System" Types

The below types are parameterized by a reference ellipse. Tags for commone ellipse's are provided, and custom ellipse's can also be used. Tag types for coordinate reference systems can be also used (e.g. OSGB36) as a parameter, but they are only a convenience for getting the coordinate reference system's ellipsoid (OSGB36 -> [AIRY ellipsoid](https://en.wikipedia.org/wiki/Ordnance_Survey_National_Grid#General))

1. `LLA`   - Longitude latitude height coordinate

2. `LL`    - Longitude latitude point coordinate to be on the ellipse's surface

3. `ECEF`  - X Y Z coordinate.  Note that ECEF has a [formal definition](https://en.wikipedia.org/wiki/ECEF) (e.g. the origin is the Earth's center of mass) which is not gauranteed to be met when transforming to this type.  More information on this provided below


The below type is parameterized by an `LL` point which determines the origin and axis directions of the local coordinate system.  The `LL` parameter can be omiited if desired.

4. `ENU`   - East North Up poistion in the local coordinate frame


### Transformation / conversion Notation

1. `convert` is used for value preserving transformations, e.g. convert(CRS{SRID(LLA{WGS84})}, LLA{WGS84}(0.0,0.0,0.0)) is appropriate as the point value remains the same.  `convert` will (should) not change values!

1. `transform` is used for value modifying transformations, e.g. transform(ECEF, LLA{WGS84}(0.0,0.0,0.0)) = ECEF{WGS84}(6.378137e6,0.0,0.0)


### Examples

##### Convert bewteen "Coordinate System" types

```julia
lla_wgs84 = LLA{WGS84}(.1167, 51.5, 0.0)

ecef_wgs84 = transform(ECEF, lla_wgs84)   # , or
ECEF(lla_wgs84)

lla2_wgs84 = LLA{WGS84}(.1168, 51.5, 0.0)
enu_ref = transform(ENU{lla_wgs84}, lla2_wgs84)    # = ENU{LLA{WGS84}(0.1167,51.5,0.0)}(6.944051663969915,4.742430141962267e-6,-3.772299267203877e-6)
enu_noref = transform(ENU, lla2_wgs84, lla_wgs84)) # = ENU{???}(6.944051663969915,4.742430141962267e-6,-3.772299267203877e-6)

lla = transform(LLA, enu_ref)
lla = transform(LLA, enu_ref, lla_wgs84)
```



##### Convert a point in a specified SRID to a "Coordinate System" Type

`srid = SRID{:EPSG, 32755}    # WGS84 datum UTM zone 55 South`
`utm = CRS{srid}(573105.43200000003, 086900.3839999996, 277.42700000000002)`
`lla_wgs84 = transform(LLA{WGS84}, utm)  # the SRID corresponding to LLA{WGS84} is known to Geodesy (see known_srids.jl).  If its not known, `
`lla_wgs84 = convert(LLA{WGS84}, transform(SRID{:EPSG, 4326}, utm))  # EPSG4326 is SRID for for the WGS84 LLA coordinate reference system`



