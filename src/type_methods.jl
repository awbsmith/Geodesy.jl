import Base.getindex


# can I make other packages add to this list?
GeodesyTypes = [LLA,
                ECEF,
                LL,
                ENU]

###########################################
# defines methods to add to point types
# Note the "geotransform" function handles all transformations
# while the "convert" function converts types leaving values unchanged
###########################################


# add some general versions of functions to overload
get_refloc{T}(::Type{T}) = error("Can't determine the reference location for type $(T)")
get_datum{T}(::Type{T}) = error("Can't determine the datum for type $(T)")
get_srid{T}(::Type{T}) = error("Can't determine the SRID for type $(T)")
get_geoid{T}(::Type{T}) = error("Can't determine the geoid for type $(T)")


# generic version of add parameters
add_param{T}(::Type{T}) = T
add_param{T}(::Type{T}, X) = T

# Hack to try to mimick v0.4's LLA.parameters style syntax
@noinline function _get_params(T::ANY, params=Type[])
    push!(params, T.var.ub)
    if isa(T.body, UnionAll)
        _get_params(T.body, params)
    end
    return params
end

###################################################################
# Function to add an integer indexing method for the geodesy types
###################################################################


# and get index methods to each type
function add_indexing{geodesy_type}(::Type{geodesy_type})
    q = quote
        getindex{T <: $(geodesy_type)}(x::T, i::Integer) = getfield(x, i)
    end
end



##############################################################
# Function to add accessor methods to Geodesy point types
##############################################################

# field input should be the field names of whatever best resemble X, Y, and Z (if it exists)
function add_accessors{geodesy_type}(::Type{geodesy_type}, fields)

    #
    # Get horizontal locations
    #
    # N.B. fields should ordered in x / y / z fashion
    q = quote

        # add horizontal accessors
        @inline getX(X::$(geodesy_type)) = X.$(fields[1])
        @inline getY(X::$(geodesy_type)) = X.$(fields[2])

        @inline get_lon(X::$(geodesy_type)) = X.$(fields[1])
        @inline get_lat(X::$(geodesy_type)) = X.$(fields[2])

        @inline get_east(X::$(geodesy_type)) = X.$(fields[1])
        @inline get_north(X::$(geodesy_type)) = X.$(fields[2])

    end

    #
    # add vertical as well
    #
    if ((length(fields)) > 2)
        qn = quote
            @inline getZ(X::$(geodesy_type)) = X.$(fields[3])
            @inline get_alt(X::$(geodesy_type)) = X.$(fields[3])
            @inline get_up(X::$(geodesy_type)) = X.$(fields[3])
        end
    else
        # for legacy support return 0 when we try to get the height of a surface point
        qn = quote
            @inline getZ(X::$(geodesy_type)) = 0.0
            @inline get_alt(X::$(geodesy_type)) = 0.0
            @inline get_up(X::$(geodesy_type)) = 0.0
        end
    end
    append!(q.args, qn.args)


    #
    # get the default parameters for this type
    #
    defaults = default_params(geodesy_type)


    #
    # get the reference location for the type if it has one
    #
    if has_refloc(geodesy_type) == Val{true}
        if (length(_get_params(geodesy_type)) == 1)  # assume its in here directly.  This will cause problems in the future
            qn = quote

                # get the reference location (note: for the moment it should always be a world surface point, but its not templated that way for convenience)
                get_refloc(::Type{$(geodesy_type)})                         =  $(defaults[1])      # it was not supplied
                get_refloc{T}(::Type{$(geodesy_type){T}})                   =  T                   # it was already supplied, but force it to be the LL type

                # get it from the point as well
                get_refloc{T}(::$(geodesy_type){T})                         =  T

            end
        else
            # need to manually define the Geodesy.get_refloc functions for this type
            qn = quote; end
        end
        append!(q.args, qn.args)
    end


    #
    # get the datum / ellipse for the type if it has one
    #
    if has_ellipse(geodesy_type) == Val{true}

        if (length(_get_params(geodesy_type)) == 1)  # it must be in here
            qn = quote

                # extract from the type
                get_datum(::Type{$(geodesy_type)})       =  UnknownDatum        # the reference point wasn't supplied so we can't get it
                get_datum{T}(::Type{$(geodesy_type){T}}) =  get_datum(T)        # it was already supplied

                # extract from the point
                get_datum{T}(::$(geodesy_type){T})       =  get_datum(T)
            end
        elseif (_get_params(geodesy_type)[1] <: AbstractDatum)  # points parameterised by datums

            qn = quote

                # extract from the type
                get_datum(::Type{$(geodesy_type)})       =  $(defaults[1])          # it was not supplied
                get_datum{T}(::Type{$(geodesy_type){T}}) =  T                       # it was already supplied

                # extract from the point
                get_datum{T}(::$(geodesy_type){T})       =  T

            end

        else
            # need to manually define the Geodesy.get_datum() functions for this type
            qn = quote; end
        end
        append!(q.args, qn.args)

        # get the actual ellipse corresponding to the above
        qn = quote
            ellipsoid{T <: $(geodesy_type)}(::Type{T})  = ellipsoid(get_datum(T))    # reference ellipsoid for the position
            ellipsoid{T <: $(geodesy_type)}(X::T)       = ellipsoid(get_datum(T))    # reference ellipsoid for the position
        end
        append!(q.args, qn.args)
    end


    #
    # get the srid for the type if it has one
    #
    if has_srid(geodesy_type) == Val{true}
        if (_get_params(geodesy_type)[1] <: AbstractSRID)
            qn = quote

                # get the parameterising datum / ellipse
                get_srid(::Type{$(geodesy_type)})       =  $(defaults[1])        # it was not supplied
                get_srid{T}(::Type{$(geodesy_type){T}}) =  T                     # it was already supplied

                # extract from the point
                get_srid{T <: $(geodesy_type)}(::T)     = get_srid(T)

                # also add the SRID contructor from this type
                (::Type{SRID})(::Type{$(geodesy_type)})         =  $(defaults[1])
                (::Type{SRID}){T <: $(geodesy_type)}(::Type{T}) =  get_srid(T)
                (::Type{SRID}){T <: $(geodesy_type)}(::T)       =  get_srid(T)
            end

            # add extra parameter stuff if this type has an extra parameter
            if (length(_get_params(geodesy_type)) > 1)
                qa = quote

                    # get the parameterising datum / ellipse
                    get_srid{T, G}(::Type{$(geodesy_type){T,G}}) =  T                     # it was already supplied

                end
                append!(qn.args, qa.args)
            end
        else
            # need to manually define the Geodesy.get_srid() function for this type
            qn = quote; end
        end
        append!(q.args, qn.args)
    end

    #
    # get the geoid for the type if it has one
    #
    if has_geoid(geodesy_type) == Val{true}
        if (_get_params(geodesy_type)[1] <: AbstractGeoid)
            qn = quote

                # get the parameterising geoid
                get_geoid(::Type{$(geodesy_type)})       =  $(defaults[1])   # it was not supplied
                get_geoid{T}(::Type{$(geodesy_type){T}}) =  T                # it was already supplied

                # extract from the point
                get_geoid{T}(::$(geodesy_type){T})       = T

            end
        elseif (_get_params(geodesy_type)[2] <: AbstractGeoid)
            qn = quote

                # get the parameterising geoid
                get_geoid(::Type{$(geodesy_type)})           =  $(defaults[2])      # it was not supplied
                get_geoid{T}(::Type{$(geodesy_type){T}})     =  $(defaults[2])      # it was not supplied
                get_geoid{T,G}(::Type{$(geodesy_type){T,G}}) =  G                   # it was supplied

                # extract from the point
                get_geoid{T,G}(::$(geodesy_type){T,G})       = G
            end
        else
            # need to manually define the Geodesy.get_geoid() function for this type
            qn = quote; end
        end
        append!(q.args, qn.args)
    end

    return q
end



##############################################################
# Function to add the template parameter to a type
##############################################################

# general version
function add_param_methods{geodesy_type}(::Type{geodesy_type})

    # get the default parameters for this type
    defaults = default_params(geodesy_type)

    # and build the code
    if (length(_get_params(geodesy_type)) == 1)
        quote
            add_param(::Type{$(geodesy_type)})        = $(geodesy_type){$(defaults[1])}    # it was not supplied
            add_param{T}(::Type{$(geodesy_type){T}})  = $(geodesy_type){T}                 # it was already supplied
        end
    elseif (length(_get_params(geodesy_type)) == 2)
        quote
            add_param(::Type{$(geodesy_type)})           = $(geodesy_type){$(defaults[1]), $(defaults[2])}    # it was not supplied
            add_param{T}(::Type{$(geodesy_type){T}})     = $(geodesy_type){T, $(defaults[2])}                 # it was half supplied
            add_param{T,U}(::Type{$(geodesy_type){T,U}}) = $(geodesy_type){T, U}                              # it was already supplied
        end
    else
        # need to manually define the Geodesy.add_param() function for this type
        quote; end
    end
end




#####################################################################
### Methods to infer template parameters based on other arguments ###
#####################################################################

#
# Function to add the template parameter to a type
#
function add_infer_param_methods{type_1, type_2}(::Type{type_1}, ::Type{type_2})

    # by default, insert the default value.  Weird.
    insert_param = Vector{Any}([x for x in default_params(type_1)])

    # and see what we can find
    for i = 1:length(_get_params(type_1))

        # what's the parameter, and does the other input have it?
        if (_get_params(type_1)[i] <: AbstractSRID) && (has_srid(type_1) == Val{true})

            # can infer the SRID
            insert_param[i] = :(get_srid(X))

        elseif (_get_params(type_1)[i] <: AbstractGeoid) && (has_geoid(type_2) == Val{true})

            # can infer the Geoid
            insert_param[i] = :(get_geoid(X))

        elseif (_get_params(type_1)[i] <: WorldPosition) && (has_refloc(type_2) == Val{true})

            # can infer the reference location
            insert_param[i] = :(get_geoid(X))

        elseif (_get_params(type_1)[i] <: AbstractDatum) && (has_ellipse(type_2) == Val{true})

            # can infer the datum
            insert_param[i] = :(get_datum(X))
        end
    end

    # and build the code
    if (length(_get_params(type_1)) == 1)
        quote
            add_param{X <: $(type_2)}(::Type{$(type_1)}, ::Type{X})            = $(type_1){$(insert_param[1])}    # it was not supplied
            add_param{T, X <: $(type_2)}(::Type{$(type_1){T}}, ::Type{X})      = $(type_1){T}                     # it was already supplied
            add_param{T <: $(type_1), X <: $(type_2)}(::Type{T}, ::X)          = add_param(T, X)
        end
    elseif (length(_get_params(type_1)) == 2)
        quote
            add_param{X <: $(type_2)}(::Type{$(type_1)}, ::Type{X})            = $(type_1){$(insert_param[1]), $(insert_param[2])}    # it was not supplied
            add_param{T, X <: $(type_2)}(::Type{$(type_1){T}}, ::Type{X})      = $(type_1){T, $(insert_param[2])}                     # it was half supplied
            add_param{T,U, X <: $(type_2)}(::Type{$(type_1){T, U}}, ::Type{X}) = $(type_1){T, U}                                      # it was already supplied
            add_param{T <: $(type_1), X <: $(type_2)}(::Type{T}, ::X)          = add_param(T, X)
        end
    else
        # need to overload Geodesy.add_param() manually
        #error("Haven't defined an add parameters method for more than 2 parameters (type: $(type_1))")
    end
end


#####################################################
# Make the constructors force typing
#####################################################

function add_constructors{geodesy_type}(::Type{geodesy_type}, fields)

    nfields = length(fields)

    q = quote

        # copy constructor
        (::Type{T}){T <: $(geodesy_type)}(X::T) = X
        (::Type{$(geodesy_type)}){T <: $(geodesy_type)}(X::T) = X
        (::Type{T}){T <: $(geodesy_type), U <: $(geodesy_type)}(X::U) = geotransform(T, X, get_handler($(geodesy_type)), get_handler($(geodesy_type)))

    end

    #
    # allow construction to convert from no datum to datum
    #

    # an expression to form a tuple from the type field
    rhs_expr = :(())
    append!(rhs_expr.args, [:(X.$(field)) for field in fields])

    # default parameters
    defaults = Geodesy.default_params(geodesy_type)
    abstracts = [supertype(def) for def in defaults]

    if (length(_get_params(geodesy_type)) == 1)
        qn = quote
            convert(::Type{$(geodesy_type){$(defaults[1])}}, X::$(geodesy_type){$(defaults[1])}) = X  # need this because it can't be resolved below
            convert{T <: $(abstracts[1])}(::Type{$(geodesy_type){T}}, X::$(geodesy_type){$(defaults[1])}) = $(geodesy_type){T}($(rhs_expr.args...))  # convert from the null parameter
            convert{T <: $(abstracts[1])}(::Type{$(geodesy_type){$(defaults[1])}}, X::$(geodesy_type){T}) = $(geodesy_type){$(defaults[1])}($(rhs_expr.args...))  # convert to the null parameter
        end
        append!(q.args, qn.args)
    elseif (length(_get_params(geodesy_type)) == 2)
        qn = quote

            # define some specifics to prevent ambiguities below
            convert(::Type{$(geodesy_type){$(defaults[1]), $(defaults[2])}}, X::$(geodesy_type){$(defaults[1]), $(defaults[2])}) = X
            convert{T <: $(abstracts[1])}(::Type{$(geodesy_type){T, $(defaults[2])}}, X::$(geodesy_type){T, $(defaults[2])}) = X
            convert{U <: $(abstracts[2])}(::Type{$(geodesy_type){$(defaults[1]), U}}, X::$(geodesy_type){$(defaults[1]), U}) = X

            # conversion to / from fully null parameters
            convert{T <: $(abstracts[1]), U <: $(abstracts[2])}(::Type{$(geodesy_type){T, U}}, X::$(geodesy_type){$(defaults[1]), $(defaults[2])}) =
                $(geodesy_type){T, U}($(rhs_expr.args...))                            # convert from the null parameter
            convert{T <: $(abstracts[1]), U <: $(abstracts[2])}(::Type{$(geodesy_type){$(defaults[1]), $(defaults[2])}}, X::$(geodesy_type){T,U}) =
                $(geodesy_type){$(defaults[1]), $(defaults[2])}($(rhs_expr.args...))  # convert to the null parameter

            # conversion to / from when the 1st parameter is known
            convert{T <: $(abstracts[1]), U <: $(abstracts[2])}(::Type{$(geodesy_type){T, U}}, X::$(geodesy_type){$(defaults[1]), U}) =
                $(geodesy_type){T, U}($(rhs_expr.args...))               # convert from the null parameter
            convert{T <: $(abstracts[1]), U <: $(abstracts[2])}(::Type{$(geodesy_type){$(defaults[1]), U}}, X::$(geodesy_type){T,U}) =
                $(geodesy_type){$(defaults[1]), U}($(rhs_expr.args...))  # convert to the null

            # conversion to / from when the 2nd parameter is known
            convert{T <: $(abstracts[1]), U <: $(abstracts[2])}(::Type{$(geodesy_type){T, U}}, X::$(geodesy_type){T, $(defaults[2])}) =
                $(geodesy_type){T, U}($(rhs_expr.args...))  # convert from the null parameter
            convert{T <: $(abstracts[1]), U <: $(abstracts[2])}(::Type{$(geodesy_type){T, $(defaults[2])}}, X::$(geodesy_type){T,U}) =
                $(geodesy_type){T, $(defaults[2])}($(rhs_expr.args...))  # convert to the null

        end
        append!(q.args, qn.args)
    end

    #
    # construction from scalars
    #


    # create an expression for a tuple x1, x2, x3...
    rhs_expr = :(())
    append!(rhs_expr.args, [:($(fields[i])) for i in 1:nfields])

    # create an expression for a tuple x1::Int, x2::Int, x3::Int... for allowing construction from Ints
    lhs_expr = :(())
    append!(lhs_expr.args, [:($(fields[i])::Real) for i in 1:nfields])

    # and build the constructor
    qn = quote
        #@compat call(::Type{$(geodesy_type)}, $(lhs_expr.args...)) = add_param($(geodesy_type))($(rhs_expr.args...))
        @compat (::Type{$geodesy_type})($(lhs_expr.args...)) = add_param($(geodesy_type))($(rhs_expr.args...))
    end
    append!(q.args, qn.args)

    # if this is an lla type, for legacy reasons we allow construction with / without the height parameter
    if has_alt(geodesy_type) == Val{true}
        if (nfields == 3)
            qn = quote  # construct it wihtout the height input
                #@compat call{T <: $(geodesy_type)}(::Type{T}, $(lhs_expr.args[1:2]...)) = add_param(T)($(rhs_expr.args[1:2]...), 0.0)
                @compat (::Type{GeodesyType}){GeodesyType <: $(geodesy_type)}($(lhs_expr.args[1:2]...)) = add_param(GeodesyType)($(rhs_expr.args[1:2]...), 0.)
            end
            append!(q.args, qn.args)
        elseif (nfields == 2)
            qn = quote  # construct it with the height input and ignore the height
                #@compat call{T <: $(geodesy_type)}(::Type{T}, $(lhs_expr.args[1:2]...), alt::Real) = add_param(T)($(rhs_expr.args[1:2]...))
                @compat (::Type{GeodesyType}){GeodesyType <: $(geodesy_type)}($(lhs_expr.args...), alt::Real) = add_param(GeodesyType)($(rhs_expr.args[1:2]...))
            end
            append!(q.args, qn.args)
        end
    end
    return q
end



#############################################################
# Add the cross constructors (constructing from other types)
#############################################################

function add_cross_constructors{type_1, type_2}(::Type{type_1}, ::Type{type_2})

    #h1, h2 = get_handler(type_1), get_handler(type_2)
    qb = quote

        # no reference point
        # call{T <: $(type_1)}(::Type{T}, X::$(type_2)) = geotransform(T, X)           # $(h1), $(h2)
        (::Type{GeodesyType}){GeodesyType <: $(type_1)}(X::$(type_2)) = geotransform(GeodesyType, X)              # $(h1), $(h2)

        # with reference point
        #call{T <: $(type_1)}(::Type{T}, X::$(type_2), ref) = geotransform(T, X, ref) # $(h1), $(h2)
        (::Type{GeodesyType}){GeodesyType <: $(type_1)}(X::$(type_2), ref) = geotransform(GeodesyType, X, ref)    # $(h1), $(h2)

    end

    # define a value preserving transform if they have different handlers
    if (get_handler(type_1) != get_handler(type_2))
        qn = quote
            convert{T <: $(type_1)}(::Type{T}, X::$(type_2)) = add_param(T)(getX(X), getY(X), getZ(X))
        end
        append!(qb.args, qn.args)
    end
    return qb
end


#####################################################
### Add Proj4 projections for each type where known
######################################################

function add_projection{geodesy_type}(::Type{geodesy_type})
    if (has_geoid(geodesy_type) != Val{true})
        quote
            get_projection{T <: $(geodesy_type)}(X::T) = get_projection(get_srid(T))
            get_projection{T <: $(geodesy_type)}(::Type{T}) = get_projection(get_srid(T))
        end
    else
        quote
            get_projection{T <: $(geodesy_type)}(X::T) = get_projection(get_srid(T), get_geoid(T))
            get_projection{T <: $(geodesy_type)}(::Type{T}) = get_projection(get_srid(T), get_geoid(T))
        end
    end
end


########################################
### Vectors of points                ###
########################################

#
# transforms
#
function add_vector_transforms{type_1, type_2}(::Type{type_1}, ::Type{type_2})

    h1, h2 = get_handler(type_1), get_handler(type_2)
    quote
        # call the geodesy vectorized version
        geotransform{T <: $(type_1), U <: $(type_2)}(::Type{T}, X::Vector{U}) = geotransform_vector(add_param(T, U), X, $(h1), $(h2))

        # call the version with a reference type
        geotransform{T <: $(type_1), U <: $(type_2)}(::Type{T}, X::Vector{U}, ref) = geotransform_vector(add_param(T, U), X, ref, $(h1), $(h2))
    end
end


#
# construction a vector of points from a matrix
#
function add_matrix_construction{geodesy_type}(::Type{geodesy_type})

    nfields = length(fieldnames(geodesy_type))

    # grab from input matrix using an index
    row_expr = :(())  # for when each point is in a row
    append!(row_expr.args, [:(X[i, $(i)]) for i in 1:nfields])

    col_expr = :(())  # for when each point is in a column
    append!(col_expr.args, [:(X[$(i), i]) for i in 1:nfields])

    quote

        # allow construction from a matrix via convert as its value preserving
        function convert{T <: $(geodesy_type)}(::Type{Vector{T}}, X::Union{AbstractMatrix, SMatrix}; row::Bool=true)
            oT = add_param(T)
            n = (row) ? size(X,1) : size(X,2)
            Xout = Vector{oT}(n)  # cant make list comprehesion get the output type right
            if row
                for i = 1:n; Xout[i] = oT($(row_expr.args...)); end
            else
                for i = 1:n; Xout[i] = oT($(col_expr.args...)); end
            end
            return Xout
        end
    end

end



#####################################################################################
# function to allow construction from fixed size vectors
# I think this is only needed because of a bug in StaticArrays
# or maybe its a new feature
#####################################################################################

function add_import_export(geodesy_type, fields)

    def_params = default_params(geodesy_type)  # default element to use
    nfields = length(fields)

    # grab from input vectors using an index
    construct_expr = :(())
    append!(construct_expr.args, [:(X[$(i)]) for i in 1:nfields])

    export_expr = :(())
    append!(export_expr.args, [:(X.$(fields[i])) for i in 1:nfields])

    q = quote

        # construct from a fixed size array vector
        convert{T <: $(geodesy_type), U <: Real}(::Type{T}, X::SVector{$(nfields), U}) = add_param(T)($(construct_expr.args...))

        # construct from a Vector
        convert{T <: $(geodesy_type), U <: Real}(::Type{T}, X::Vector{U}) = add_param(T)($(construct_expr.args...))

        # construct from a tuple
        convert{T <: $(geodesy_type), U}(::Type{T}, X::NTuple{$(nfields), U}) = add_param(T)($(construct_expr.args...))

        # export to a fixed size array vector
        convert(::Type{SVector}, X::$(geodesy_type)) = SVector{$(nfields), Float64}($(export_expr.args...))
        convert{U <: Real}(::Type{SVector{$(nfields), U}}, X::$(geodesy_type)) = Vec{$(nfields), U}($(export_expr.args...))

        (::Type{SVector})(X::$(geodesy_type)) = convert(SVector, X)
        (::Type{SVector{$(nfields), U}}){U <: Real}(X::$(geodesy_type)) = convert(SVector{$(nfields), U}, X)

        # export to a vector
        convert(::Type{Vector}, X::$(geodesy_type)) = vcat($(export_expr.args...))

        # export to a tuple
        convert(::Type{NTuple}, X::$(geodesy_type)) = $(export_expr)

    end
    return q
end




#
# Function build a code block with everything from above
# The idea is other packages can use this to create code blocks,
# and evaluate it in the Geodesy module to save exporting everything.  IDK
#
function build_methods{geo_type}(::Type{geo_type}, fields=fieldnames(geo_type), known_types=[])

    # create a code block to amalgamate all the othe code blocks
    qb = Expr(:block)

    # add indexing method
    append!(qb.args, Geodesy.add_indexing(geo_type).args)

    # add accessor functions
    append!(qb.args, Geodesy.add_accessors(geo_type, fields).args)

    # add the parameter manipulation methods
    append!(qb.args, Geodesy.add_param_methods(geo_type).args)

    # add basic maths
    append!(qb.args, Geodesy.add_maths(geo_type, fields).args)

    # add additional constructors
    append!(qb.args, Geodesy.add_constructors(geo_type, fields).args)

    # the number of elements in the representation
    append!(qb.args, Geodesy.add_projection(geo_type).args)

    # add constructing vectors of these things from matrices
    append!(qb.args, Geodesy.add_matrix_construction(geo_type).args)

    # add creating them from Vecs, Vectors, and tuples
    append!(qb.args, Geodesy.add_import_export(geo_type, fields).args)

    # NaN checks
    # append!(qb.args, Geodesy.add_nan_check(geo_type).args)

    # and add the binary methods
    for other_type in known_types

        # build both ways
        for comb in [(geo_type, other_type), (other_type, geo_type)]

            # add the parameter inference methods (both ways around)
            append!(qb.args, Geodesy.add_infer_param_methods(comb[1], comb[2]).args)

            # add the construction from other geodesy type methods (both ways around)
            append!(qb.args, Geodesy.add_cross_constructors(comb[1], comb[2]).args)

            # add vectorized transforms between types
            append!(qb.args, Geodesy.add_vector_transforms(comb[1], comb[2]).args)

        end
    end

    # add this new point to the list of types
    push!(known_types, geo_type)

    return qb
end


########################################
### Finally, make everything!        ###
########################################


#
# Now go through and add everything
#
for (i, t) in enumerate(GeodesyTypes)
     Geodesy.eval(build_methods(t, fieldnames(t), copy(GeodesyTypes[1:i-1])))  # TODO: eval? figure out the right way to do this
end

