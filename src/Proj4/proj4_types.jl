# add a utm point type since we can handle it with Proj4
# leverage the point type construction code used in the Geodesy package
eval(BuildPointType(:UTM,  UTM_CS,  (:false_east, :false_north, :up), nothing))


# allow conversion constructing it
# add self conversions
    eval(
        quote
            @inline convert{T <: $(ptype)}(::Type{T}, X::T) = X
            @inline convert(::Type{$(ptype)}, X::$(ptype)) = X
        end
        )

    # add conversions to other types
    for otype in setdiff(subtypes(ABSTRACT_POSITION), [UnknownPoint, ptype])
        eval(
            quote
                @inline convert{Tdest <: $(otype)}(::Type{Tdest}, X::$(ptype)) = geotransform(Tdest, X)
                @inline convert{Tdest <: $(otype)}(::Type{Tdest}, X::$(ptype), crs) = geotransform(Tdest, get_crs(Tdest), X, crs)
            end
            )
    end

