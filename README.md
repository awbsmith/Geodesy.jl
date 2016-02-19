# Geodesy

[![Build Status](https://travis-ci.org/JuliaGeo/Geodesy.jl.svg?branch=master)](https://travis-ci.org/JuliaGeo/Geodesy.jl)
[![Coverage Status](http://img.shields.io/coveralls/JuliaGeo/Geodesy.jl.svg)](https://coveralls.io/r/JuliaGeo/Geodesy.jl)

Work with points defined in various coordinate systems. The code has been split out from [OpenStreetMap.jl](https://github.com/tedsteiner/OpenStreetMap.jl), and functionality expanded.

Coordinate systems `LL`, `LLA`, `ECEF`, `SRID_Pos`, `ENU` are supported. Transforms between between those types are supported.

While `LL`, `LLA`, and `ECEF` are parameterized by a reference ellipsoid (not a datum!) and are the intended type to with while using this package. 

The `SRID_Pos` type (parameterized by an `SRID` which identifies the coordinate system and datum) is provided to aid importing and exporting data from various datums. Transformations involving an `SRID\_Pos` point are performed by the [Proj4](https://github.com/FugroRoames/Proj4.jl) package.

The `ENU` local point type may be parameterized by an `LL` point describing the origin of the coordinate system if desired

Example usage is shown in example_usage.jl

The full list of types, constants, and methods provided is at the top of [src/Geodesy.jl](src/Geodesy.jl).
