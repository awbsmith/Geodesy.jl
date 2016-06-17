import Base: convert

# list of srids for known datums (allows converting LLA / ECEF to SRID points)
# TODO: Check if these are actually correct!

# defines whether we want the below dictionary available via dispath instead of the dict.
# This allows the "correct" point container to be chosen (with type safety)
# but also clogs up the type system
create_typed_srid_conversion = true

immutable SRID_Dict_Type
    SRID_to_CRS::Dict{SRID, CRS{DataType, DataType}} # N.B. Using datatype here for type safety
    CRS_to_SRID::Dict{CRS{DataType, DataType}, SRID}
end

function build_srid_dict()

    # TODO: Parse Proj4's list to build this?

    # build the list
    Dict(

    ############
    # WGS84
    ############

    SRID(:EPSG, 4326) => CRS{DataType, DataType}(LLA_CS, WGS84),
    SRID(:EPSG, 4978) => CRS{DataType, DataType}(ECEF_CS, WGS84),

    ##############
    # NAD27
    ##############

    SRID(:EPSG, 4267) => CRS{DataType, DataType}(LLA_CS, NAD27),

    ##############
    # ED50
    ##############

    SRID(:EPSG, 4230) => CRS{DataType, DataType}(LLA_CS, ED50),

    ##############
    # OSGB36
    ##############

    SRID(:EPSG, 4277) => CRS{DataType, DataType}(LLA_CS, OSGB36),

    ##############
    # GDA94
    ##############

    # Australia
    SRID(:EPSG, 4283) => CRS{DataType, DataType}(LLA_CS, GDA94),

    ##############
    # ETRS89
    ##############

    # Europia
    SRID(:EPSG, 4258) => CRS{DataType, DataType}(LLA_CS, ETRS89),

    ##############
    # NAD83
    ##############

    # Americania
    SRID(:EPSG, 4269) => CRS{DataType, DataType}(LLA_CS, NAD83)

    )
end

# now build the reverse dictionary
function build_crs_dict()
    srid_dict = build_srid_dict()
    crs_dict = Dict{CRS{DataType, DataType}, SRID}()
    for key in keys(srid_dict)
        crs_dict[srid_dict[key]] = key
    end
    return crs_dict
end
SRID_Dict = SRID_Dict_Type(build_srid_dict(), build_crs_dict())


# get the SRID for a CRS
function get_srid{CS, DATUM}(crs::CRS{CS, DATUM})
    crs_fetch = CRS{DataType, DataType}(CS, DATUM)
    if haskey(SRID_Dict_Type.SRID_to_CRS, crs_fetch)
        srid = SRID_Dict_Type.SRID_to_CRS[crs_fetch]
    else
        error("The SRID is unknown for crs: $(crs)")
    end
end
convert(::Type{SRID}, crs::CRS) = get_srid(crs)


# get the CRS of an SRID.  N.B. this returns a CRS{DataType,DataType} for type safety
function get_crs(srid::SRID)
    if haskey(SRID_Dict_Type.CRS_to_SRID, srid)
        crs = SRID_Dict_Type.CRS_to_SRID[srid]
    else
        error("The CRS is unknown for srid: $(crs)")
    end
end
convert(::Type{CRS}, srid::SRID) = UnknownCRS  # nothing else is type safe

#
# Repeat the above for the tagged SRIDs
# The switch is for whether we want to declare a new type for each SRID in the dictionary
# doing so allows more type inference, but clogs up the type system otherwise
#
if (create_typed_srid_conversion)

    #
    # create the typed versions for dispatching onto
    #
    qb = quote end
    for crs in keys(SRID_Dict.CRS_to_SRID)
        srid = SRID_Dict.CRS_to_SRID[crs]
        auth_expr = Expr(:quote, srid.auth)
        qn = quote
            @compat (::Type{SRIDt})(::CRS{$(crs.cs), $(crs.datum)}) = SRIDt{$(auth_expr), $(srid.code)}(nothing, nothing)
            @compat (::Type{CRS})(::SRIDt{$(auth_expr), $(srid.code)}) = CRS($(crs.cs), $(crs.datum))
        end
        append!(qb.args, qn.args)
    end
    eval(qb)

else

    # convert the typed SRID to a data one instead and go
    @inline convert{AUTH, CODE}(::Type{CRS}, ::SRIDt{AUTH, CODE}) = CRS(SRID(AUTH, CODE))
    @inline convert(::Type{SRIDt}, crs::CRS) = SRIDt(SRID(crs))

end

