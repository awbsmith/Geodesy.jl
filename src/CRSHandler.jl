
# define handlers for different CRS's
abstract CRS_Handler

# geodesy function
immutable GeodesyHandler <: CRS_Handler; end

# CoordinateTransforms style
immutable CoordinateTransformsHandler <: CRS_Handler; end

# default is to try and handle conversions usign this package
# T1 and T2 are expected to be CRS's
@inline get_handler{T1, T2}(::T1, ::T2) = GeodesyHandler

