# Geodesy
[![Build Status](https://travis-ci.org/JuliaGeo/Geodesy.jl.svg?branch=master)](https://travis-ci.org/JuliaGeo/Geodesy.jl)
[![Coverage Status](http://img.shields.io/coveralls/JuliaGeo/Geodesy.jl.svg)](https://coveralls.io/r/JuliaGeo/Geodesy.jl)

The Geodesy package implements geodetic transformations bewteen various coordinate systems. 

### Quick Guide

```julia
# start with a GPS coordinate, which is a [longitude, latitude, altitude] (LLA) location in the WGS84 datum
lla_src = LLA{WGS84}(.1167, 51.5, 0.0)
# = Geodesy.LLA{WGS84} 
#     lon: Float64 0.1167
#     lat: Float64 51.5
#     alt: Float64 0.0


# define a destination GPS coordinate
lla_dest = LLA{WGS84}(0.11706002390375786,51.500449405987546,0.00024498533457517624)

# convert them to Earth Centered Earth Fixed corrdinates (a Cartesian cordinate system)
ecef_src = ECEF(lla_src)

# transform the destination into a coordinate frame centered at the source location, and with axis directions East, North, and Up (ENU)
# (this basically gives a direction from lla_src to lla_dest)
enu_dest = ENU(lla_dest, lla_src)
# = Geodesy.ENU{Geodesy.UnknownRef} 
#        east:  Float64 25.00000000000066
#        north: Float64 49.99999999967102
#        up:    Float64 -3.318305630273244e-10

# show the datums and subtypes of datums known to this package
get_datums()
#=    "ED50 (StaticDatum)"                   
#     "NAD27 (StaticDatum)"                  
#     "OSGB36 (StaticDatum)"                 
#     "WGS84 (StaticDatum)"                  
#     "ETRS89 (DynDatum)"                    
#     "GDA94 (DynDatum)"                     
#     "NAD83 (DynDatum)"                     
#     "UnknownDatum"
#     "AIRY_ELLIPSE (Ellipse)"               
#     "CLARKE66_ELLIPSE (Ellipse)"           
#     "CustomEllipse{T<:Ellipsoid} (Ellipse)"
#     "GRS80_ELLIPSE (Ellipse)"              
#     "HAYFORD_ELLIPSE (Ellipse)"            
#     "WGS84_ELLIPSE (Ellipse)"              
```
More detailed usage is given below


### Terminology

![Coordinate Reference systems](http://www.crs-geo.eu/SharedDocs/Bilder/CRS/schema-crs-datum-cs,property=default.gif) 
[www.crs-geo.eu](http://www.crs-geo.eu/nn_124224/crseu/EN/CRS__Overview/definition-crs__node.html). 

[View the above as a data structure](http://i.stack.imgur.com/aeS8k.png)

The above images gives a quick picture of the components of the coordinate reference systems used in geodesy, with the data structure link showing the data structure used by the EPSG authority for the WGS84 (GPS) coordinate reference system as shown on their [EPSG website](http://www.epsg-registry.org/)).

This Geodesy package is intended for use with the "Coordinate System" subtypes shown above. Following this Geodesy position types do not have full datum knowledge (e.g. where the coordinate system's origin is relative to the Earth), only the reference ellipsoid required to perform the coordindate system transforms.  Transforms defined in this package convert between coordinate systems (e.g. longitude latitude height -> cartesian etc), although a coordinate reference system type is currently provided to facilitate importing and exporting data.


### "Coordinate System" Types

The below types are parameterised by a reference datum. Note that the coordinate system types only "understands" the datum's ellipse; however using a datum type as a parameter is a convenient way to get the reference ellipse while also stopping direct comparison of points from different datums (that may use the same ellipse).  A discussion on datums vs ellipses is given later in this readme.

Some common ellipses and datums are provided, and custom ellipse's can also be used. For a list of all predefined datums, use `Geodesy.get_datums()`

1. `LLA`   - Longitude, [(geodetic) latitude](https://en.wikipedia.org/wiki/Latitude#Geodetic_and_geocentric_latitudes), altitude coordinate

2. `LL`    - Longitude, [(geodetic) latitude](https://en.wikipedia.org/wiki/Latitude#Geodetic_and_geocentric_latitudes) coordinate to be on the ellipse's surface

3. `ECEF`  - Cartesian coordinate.  Note that ECEF has a [meaning beyong a Cartesian coordinate system](https://en.wikipedia.org/wiki/ECEF), however it is used in this package to mean only a Cartesian coordinate system with origin equal to the ellipse's center.  

The below type is parameterized by an `LL` point which determines the origin and axis directions of the local coordinate system.  The `LL` parameter can be omitted from the template and passed to function seperately if desired.

4. `ENU`   - East North Up position in the local coordinate frame.


### "Coordinate **Reference** System" Types

Coordinate reference sytem types have knowledge of the datum which is required to map a known point on the Earth to a position in the coordinate system. The [Proj4 package](https://github.com/FugroRoames/Proj4.jl) is currently used as a backend to allow transforming to / from / between "Coordinate _Reference_ Systems". 

The below type is parameterised by a [spatial reference ID (SRID)](https://en.wikipedia.org/wiki/SRID). Note that SRIDs are used to describe more than just coordinate reference systems, so take care when selecting a suitable SRID.

1. `CRS` - The coordinate _reference_ system point type.  This type should be used for operations that require knowledge of the datum, e.g. swapping between datums. Transformations involving this type are perfromed by Proj4 as a full understanding of coordinate reference systems is outside the scope of this package. It's left to the user to correctly interpret the fields of this type as different SRIDs use different coordinate systems ([lat, long, height] / [x, y, z] / [false east, false north, height] / etc).

**Roadmap note**: Coordinate reference system types should be migrated to backend packages (Proj4 etc), with the Geodesy package providing common interface methods for them (functions to overload etc).  This change would mean Proj4 depends on Geodesy instead of Geodesy depending on Proj4 matching Geodesy being a "low level" package.  This also allows different backend packages to have their own coordinate reference system type.


#### Datums vs Ellipsoids

Datums contain more information than just the reference ellipse, they also contain the information to map a fixed point on the Earth's surface to a point on the ellipse.  Two different datums can use the same ellipse, but have a different origin and / or orientation relative to a set of points on the Earth's surface because they use different mappings from the ellipse to the Earth's surface.  For example, the _WGS84 datum_ ([EPSG6326](https://epsg.io/6326-datum) as used by GPS systems) and the _Vietnam2000 datum+ ([EPSG4756](https://epsg.io/4756-5194)) both use the _WGS84 ellipsoid_ ([EPSG7030](https://epsg.io/7030-ellipsoid)) but will have a different coordinate for the same place on Earth (e.g. Ho Chi Minh City), because the ellipses have been aligned relative to the Earth in different ways.  In the Geodesy package, 

Comparing points in two different datums requires knowing the transformation to align the ellipses.  This information is only known for the "coordinate _reference_ system" (`CRS`) type.



### Transformations vs Conversions

1. `geotransform` is used for value modifying transformations, e.g.:
```julia
    lla_wgs84 = LLA{WGS84}(0.0,0.0,0.0)
    geotransform(ECEF, lla_wgs84) # = ECEF{WGS84}(6.378137e6,0.0,0.0)
```


2. `convert` is used for value preserving transformations, e.g. change the position's type to include the coordinate reference system:
```julia
    lla_wgs84 = LLA{WGS84}(0.0,0.0,0.0)  # position in the LLA coordinate system using the WGS84 reference ellipsoid.  Note multiple datums use the WGS84 ellipse.
    lla_wgs84_srid = SRID(LLA{WGS84})    # define the coordinate reference system for the WGS84 datum with LLA coordinate system

    # the position is the same!    
    convert(CRS{lla_wgs84_srid}, lla_wgs84) # = CRS{EPSG4326}(0.0,0.0,0.0)    
```
`convert` will / should not change values!

Type constructors are overloaded to use the `geotransform` function where applicable


### Examples

##### Basic Usage

```julia

# start with an [longitude, latitude, altitude] (LLA) location in the WGS84 datum (i.e. a GPS coordinate)
lla_src = LLA{WGS84}(.1167, 51.5, 0.0)

# define a destination GPS coordinate
lla_dest = LLA{WGS84}(0.11706002390375786,51.500449405987546,0.00024498533457517624)

# Transform them both into a Cartesian Earth Centered Earth Fixed (ECEF) frame
ecef_src = ECEF(lla_src)
# = Geodesy.ECEF{WGS84} 
#     x: Float64 3.9786402778542214e6
#     y: Float64 8103.702688750244
#     z: Float64 4.968362457291028e6
ecef_dest = ECEF(lla_src)


# transform the destination into a coordinate frame centered at the source location, and with axis direction East, North, and Up (ENU)
# (this gives a direction from lla_src to lla_dest)
enu_dest = ENU(lla_dest, lla_src)
# = Geodesy.ENU{Geodesy.UnknownRef} 
#        east:  Float64 25.00000000000066
#        north: Float64 49.99999999967102
#        up:    Float64 -3.318305630273244e-10

```

##### Advanced Usage 


###### Use a "Coordinate Reference System" type to represent a coordinate reference system specified by an SRID

```julia

# create a point in a coordinate reference system specified by an SRID
srid = SRID{:EPSG, 32755}                                                        # [WGS84 datum UTM zone 55 South](https://epsg.io/32755)
X_utm = CRS{srid}(573105.43200000003, 086900.3839999996, 277.42700000000002)     # a point in the above coordinate system

# convert to a "known" type
lla_wgs84 = geotransform(LLA{WGS84}, X_utm)                                      # the SRID corresponding to LLA{WGS84} is known to Geodesy (see known_srids.jl).  

# the above line is a convenient way to perform the steps
Xo = geotransform(CRS{SRID{:EPSG, 4326}}, X_utm)                                 # SRID{:EPSG, 4326} corresponds to [LLA with a WGS84 datum](https://epsg.io/4326)
lla_wgs84 = convert(LLA{WGS84}, Xo)                                              # a value presere type conversion


```

###### Change the backend package used to perform a conversion

```julia

    # do a geotransformation the usual way
    lla = LLA{WGS84}(.1167, 51.5, 0.0)
    ecef = geotransform(ECEF{WGS84}, lla)

    # now make the Proj4 package do the conversion, by specifying the type handlers in the function inputs 
    # Note that Proj4 will perfrom the transform is either src or dest used the Proj4 Handler
    ecef_proj4 = geotransform(ECEF{WGS84}, lla, Geodesy.Proj4Handler, Geodesy.Proj4Handler)  # order is dest handler, then source handler

    # how different is the result?
    dX = ecef - ecef_proj4
    @printf("%0.19g", norm(dX))

```


###### Exploring Template Parameter Usage

The usual style in this package is to use template parameters to supply an extra argument to the `geotransform` function.  

Template parameters for the different types are useful in that they save passing information throughout your code, but can cause issues surrounding type stability, or when you use so many different template parameters that you'll kill Julia's type system.

####### Datum template parameter example:

```julia

    # create an LLA points with a known datum
    lla_known_datum = LLA{WGS84}(.1167, 51.5, 0.0)

    # the transformation to ECEF require parameters for the datum's ellipse
    # When we use the templated version: 
    ecef_known_datum = geotransform(ECEF, lla_known_datum)
    # Geodesy.ECEF{WGS84} 
    #      x: Float64 3.9786402778542214e6
    #      y: Float64 8103.702688750244
    #      z: Float64 4.968362457291028e6
    #
    # the ellipse parameters are added to the function inputs, making the above line equivilent to:
    ecef_known_datum = geotransform(ECEF, lla_known_datum, Geodesy.ellipsoid(typeof(lla_known_datum)))

    # we can use non-templated version by specifying the ellipse manually
    lla_unknown_datum = LLA(.1167, 51.5, 0.0)
    ellipse = Geodesy.eWGS84  # WGS84 ellipse parameters (N.B. Geodesy.Ellipsoid(...) can be used to create custom ellipsoids)
    ecef_unknown_datum = geotransform(ECEF, lla_unknown_datum, ellipse)
    # = Geodesy.ECEF{Geodesy.UnknownDatum} 
    #      x: Float64 3.9786402778542214e6
    #      y: Float64 8103.702688750244
    #      z: Float64 4.968362457291028e6


    # and check the result
    using FixedSizeArrays
    dX = Vec(ecef_known_datum) - Vec(ecef_unknown_datum)  # convert to Vec because they're different types
    @printf("%0.19g", norm(dX))

```

####### Reference point template parameter example:

```julia

    # create two LLA points
    lla_ref = LLA{WGS84}(.1167, 51.5, 0.0)
    lla_dest = LLA{WGS84}(0.11706002390375786,51.500449405987546,0.00024498533457517624)

    # when we transform to a local reference frame, the default behaviour is not to add the reference point to the output type
    # so the reference point must be added to the function call
    enu_noref = geotransform(ENU, lla_dest, lla_ref)
    # = Geodesy.ENU{Geodesy.UnknownRef} 
    #      east: Float64 25.00000000000066
    #      north: Float64 49.99999999967102
    #      up: Float64 -3.318305630273244e-10

    
    # we can specify the reference point in the ENU type
    enu_ref = geotransform(ENU{lla_ref}, lla_dest)
    # = Geodesy.ENU{Geodesy.LLA{WGS84}(0.1167,51.5,0.0)} 
    #             east: Float64 25.00000000000066
    #             north: Float64 49.99999999967102
    #             up: Float64 -3.318305630273244e-10

    # which means we don't have to supply the reference point in later conversions
    lla_out = geotransform(LLA, enu_ref)

```


####### SRID template parameter example:

```julia

    # start off by defining a source and destination coordinate reference system by their SRID
    srid_src = SRID{:EPSG, 4326}   # corresponds to LLA{WGS84}
    srid_dest = SRID{:EPSG, 4978}  # corresponds to ECEF{WGS84}

    # and create a point
    X = CRS{srid_src}(.1167, 51.5, 0.0)

    # when we perform the transformation, the SRID information gets added as an extra parameter to the transformtion function,
    # so
    Xo = geotransform(CRS{srid_dest}, X)
    # = Geodesy.CRS{EPSG4978} 
    #    x: Float64 3.9786402778542214e6
    #    y: Float64 8103.702688750244
    #    z: Float64 4.968362457291028e6
    #
    # is equivilent to: 
    Xo = geotransform(CRS{srid_dest}, X, (srid_dest, srid_src))
    # = Geodesy.CRS{Geodesy.UnknownSRID} 
    #     x: Float64 3.9786402778542214e6
    #     y: Float64 8103.702688750244
    #     z: Float64 4.968362457291028e6

    
    # so we can do the above without using the template parameters
    X_noref = CRS(.1167, 51.5, 0.0)
    Xo_noref = geotransform(CRS, X_noref, (srid_dest, srid_src))
    # = Geodesy.CRS{Geodesy.UnknownSRID} 
    #     x: Float64 3.9786402778542214e6
    #     y: Float64 8103.702688750244
    #     z: Float64 4.968362457291028e6

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
    custom_ecef = geotransform(CustomECEF, custom_lla)      # = CustomECEF(2016-02-25T16:13:59,"test input",3.9786402778542214e6,8103.702688750244,4.968362457291028e6)

```




####### Datums vs Ellipses Example 

```julia

# The [vn2000 datum](https://epsg.io/4756-5194) datum for Vietnam uses the WGS84 ellipsoid 
vn2000 = SRID{:EPSG, 4756}  # SRID to parameterize the Coordinate Reference System

# define this datum for use with Geodesy
immutable VN2000 <: Geodesy.Datum; end
import Geodesy: ellipsoid
ellipsoid(::Type{VN2000}) = ellipsoid(Geodesy.WGS84_ELLIPSE)  # overload the function to retrieve this datum's ellipsoid

# define an LLA point in this datum
lla_vn2000 = CRS{vn2000}(0, 0, 0)

# convert this point to the WGS84 datum
lla_wgs84 = geotransform(LLA{WGS84}, crs_vn2000)              # = LLA{WGS84}(-0.0003528546332018209,-0.0010055677147903898,-192.75199438724667)

```

We can check the Proj4 well known text for the VN2000 datum to see the *+towgs84* transformation parameters

```Julia
Geodesy.proj4_str(vn2000) # = "+proj=longlat +ellps=WGS84 +towgs84=-192.873,-39.382,-111.202,-0.00205,-0.0005,0.00335,0.0188 +no_defs"
```


#### A note on the ECEF type

An ECEF coordinate system's origin should be the Earth's center of mass and have axes aligned with the International Reference Pole and International Reference Meridian ([ECEF](https://en.wikipedia.org/wiki/ECEF)).  Since the coordinate system types only use the ellipse part of the datum, they have no information to align the datum's axes to the International Reference Meridian etc. 

As a result the ECEF point type used in this package is simply a Cartesian coordinate system with an origin set to the parameterising ellipse's center, with axes specified by the ellipse which could point anywhere on the Earth's surface.  The tag ECEF is choosen for pragmatism.

As an example:

```julia
lla_osgb36 = LLA{Geodesy.OSGB36}(0, 0, 0)      # OSGB36 is a good match to the Earth in the UK but not elsewhere, and is not centered on the Earth's center of mass (i.e. not a [geocentric datum](http://support.esri.com/en/knowledgebase/GISDictionary/term/geocentric%20datum))
ecef_geodesy = geotransform(ECEF, lla_osgb36)
# = Geodesy.ECEF{OSGB36} 
#     x: Float64 6.377563396e6
#     y: Float64 0.0
#     z: Float64 0.0
                      

# we can compare the above to the ECEF position as calculated by the Proj4 package
ecef_proj4 = geotransform(ECEF{WGS84}, lla_osgb36, Geodesy.Proj4Handler, Geodesy.Proj4Handler)                      # convert to a type native to this package (type conversion not a transformation)
# = Geodesy.ECEF{WGS84} 
#      x: Float64 6.377879171552554e6
#      y: Float64 -99.12039106890559
#      z: Float64 534.423089412207


# the difference
dist = norm(Vector(cart) - Vector(ecef))                   # = 628.6072621385788



```














