import Base.show

strip_type_params{pType}(::Type{pType}) = string(pType.name)
strip_module_name{pType}(::Type{pType}) = strip_module_name(string(pType))
strip_module_name(str::AbstractString) = replace(str, "^.*\.", "")

function minimize_display{pType}(::Type{pType})

    # get rid of the module identifier
    # raw_str = strip_module_name(strip_type_params(pType))
    raw_str = replace(string(pType), "Geodesy2.", "")
    # raw_str = strip_module_name(pType)
    raw_str = replace(raw_str, r"{.*}", "")
    raw_expr = Symbol(raw_str)

    # create a data format string
    fmt_str = "(%s"
    for i in 2:ndims(pType)
        fmt_str *= ", %s"
    end
    data_args = [:(string(X[$i])) for i in 1:ndims(pType)]

    # build the quote block
    quote

        # show the type
        show(io::Base.IO, ::Type{$(raw_expr)}) = print(io, $(raw_str))

        # only the CRS in the type
        function show{CS, DATUM}(io::Base.IO, ::Type{$(raw_expr){CRS{CS, DATUM}}})
            if (CS != $(default_coord_system(pType))) # what's the user up to?
                print(io, $(raw_str) * "{" * string(CRS{CS, DATUM}) * "}")
            elseif (DATUM == UnknownDatum)
                print(io, $(raw_str))
            else
                print(io, $(raw_str) * "{" * string(DATUM) * "}")
            end
        end

        # CRS and element type in the type
        function show{CS, DATUM, eT}(io::Base.IO, ::Type{$(raw_expr){CRS{CS, DATUM}, eT}})
            if (CS != $(default_coord_system(pType))) # what's the user up to?
                print(io, $(raw_str) * "{" * string(CRS{CS, DATUM}) * ", " * string(eT) * "}")
            elseif DATUM == UnknownDatum
                print(io, $(raw_str) * "{" * string(eT) * "}")
            else
                print(io, $(raw_str) * "{" * string(DATUM) * ", " * string(eT) * "}")
            end
        end

        # display the data as well
        function show{CRS, eT}(io::Base.IO, X::$(raw_expr){CRS, eT})
            if (sizeof(CRS) > 0) || (get_coord_system(X.crs) != $(default_coord_system(pType)()))
                print(io, string($(raw_expr){CRS, eT}) * @sprintf($(fmt_str * ", %s)"),  $(data_args...), string(X.crs)))  # it can't be reconstructed from the type so show the crs
            else
                print(io, string($(raw_expr){CRS, eT}) * @sprintf($(fmt_str * ")"),  $(data_args...)))
            end
        end
    end
end

# codegen display methods
for pType in [LLA, LL, ECEF, ENU]
    eval(minimize_display(pType))
end


#
# Remove the module name when displaying coordinate systemss
#

function minimize_cs_display{csType}(::Type{csType})
    raw_str = replace(string(csType), "Geodesy2.", "")
    qb = quote
        show(io::Base.IO, ::Type{$(csType)}) = print(io, $(raw_str))
        show(io::Base.IO, ::$(csType)) = print(io, $(raw_str) * "()")
    end
end


for csType in get_CS_types()
    eval(minimize_cs_display(csType))
end

