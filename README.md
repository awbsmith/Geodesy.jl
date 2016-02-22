# Geodesy

#[![Build Status](https://travis-ci.org/JuliaGeo/Geodesy.jl.svg?branch=master)](https://travis-ci.org/JuliaGeo/Geodesy.jl)
#[![Coverage Status](http://img.shields.io/coveralls/JuliaGeo/Geodesy.jl.svg)](https://coveralls.io/r/JuliaGeo/Geodesy.jl)

Work with points defined in various coordinate systems. The code has been split out from [OpenStreetMap.jl](https://github.com/tedsteiner/OpenStreetMap.jl), and functionality expanded.

Coordinate systems `LL`, `LLA`, `ECEF`, `SRID_Pos`, `ENU` are supported. Transforms between between those types are supported.

While `LL`, `LLA`, and `ECEF` are parameterized by a reference ellipsoid (not a datum) and are the intended type to work with while using this package. 

The `SRID_Pos` type (parameterized by an `SRID` which identifies the coordinate system and datum). Transformations involving an `SRID_Pos` point are performed by the [Proj4](https://github.com/FugroRoames/Proj4.jl) package.

The `ENU` local point type may be parameterized by an `LL` point describing the origin of the coordinate system if desired

The intended workflow when using this package is to:
1) import data using the `SRID_Pos` point type (if not already in a native type)
2) transform these points into this package's native types (`LL`, `LLA`, `ECEF`, `ENU`),
3) Process the data
4) export data ausing `SRID_Pos` point type (if a non native export type is required)

Example usage is shown in example_usage.jl


