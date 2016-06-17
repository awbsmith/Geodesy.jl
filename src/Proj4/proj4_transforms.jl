


# use Proj4 as the default handler for points specified by a CRS
get_handler{T <: }(::Type{T}) = Proj4Handler
