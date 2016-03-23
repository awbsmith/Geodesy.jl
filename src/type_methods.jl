GeodesyTypes = [LLA,
                ECEF,
                CRS,
                CCRS_Geoid,
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
        if (length(geodesy_type.parameters) == 1)  # assume its in here directly.  This will cause problems in the future
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

        if (geodesy_type.parameters[1] <: AbstractDatum)  # points parameterised by datums 
            qn = quote
            
                # extract from the type
                get_datum(::Type{$(geodesy_type)})       =  $(defaults[1])          # it was not supplied
                get_datum{T}(::Type{$(geodesy_type){T}}) =  T                       # it was already supplied

                # extract from the point
                get_datum{T}(::$(geodesy_type){T})       =  T

            end
        elseif (length(geodesy_type.parameters) == 1)  # it must be in here
            qn = quote

                # extract from the type
                get_datum(::Type{$(geodesy_type)})       =  UnknownDatum        # the reference point wasn't supplied so we can't get it
                get_datum{T}(::Type{$(geodesy_type){T}}) =  get_datum(T)        # it was already supplied

                # extract from the point
                get_datum{T}(::$(geodesy_type){T})       =  get_datum(T) 
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
        if (geodesy_type.parameters[1] <: AbstractSRID)  
            qn = quote
                
                # get the parameterising datum / ellipse
                get_srid(::Type{$(geodesy_type)})       =  $(defaults[1])        # it was not supplied
                get_srid{T}(::Type{$(geodesy_type){T}}) =  T                     # it was already supplied

                # extract from the point
                get_srid{T <: $(geodesy_type)}(::T)     = get_srid(T)

                # also add the SRID contructor from this type
                call(::Type{SRID}, ::Type{$(geodesy_type)})         =  $(defaults[1])
                call{T <: $(geodesy_type)}(::Type{SRID}, ::Type{T}) =  get_srid(T)
                call{T <: $(geodesy_type)}(::Type{SRID}, ::T)       =  get_srid(T)
            end
            
            # add extra parameter stuff if this type has an extra parameter
            if (length(geodesy_type.parameters) > 1)
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
        if (geodesy_type.parameters[1] <: AbstractGeoid)  
            qn = quote

                # get the parameterising geoid
                get_geoid(::Type{$(geodesy_type)})       =  $(defaults[1])   # it was not supplied
                get_geoid{T}(::Type{$(geodesy_type){T}}) =  T                # it was already supplied

                # extract from the point
                get_geoid{T}(::$(geodesy_type){T})       = T

            end
        elseif (geodesy_type.parameters[2] <: AbstractGeoid)
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
    if (length(geodesy_type.parameters) == 1)
        quote
            add_param(::Type{$(geodesy_type)})        = $(geodesy_type){$(defaults[1])}    # it was not supplied
            add_param{T}(::Type{$(geodesy_type){T}})  = $(geodesy_type){T}                 # it was already supplied
        end
    elseif (length(geodesy_type.parameters) == 2)
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
    for i = 1:length(type_1.parameters)

        # what's the parameter, and does the other input have it?
        if (type_1.parameters[i] <: AbstractSRID) && (has_srid(type_1) == Val{true})

            # can infer the SRID
            insert_param[i] = :(get_srid(X))

        elseif (type_1.parameters[i] <: AbstractGeoid) && (has_geoid(type_2) == Val{true})

            # can infer the Geoid    
            insert_param[i] = :(get_geoid(X))

        elseif (type_1.parameters[i] <: WorldPosition) && (has_refloc(type_2) == Val{true})

            # can infer the reference location
            insert_param[i] = :(get_geoid(X))

        elseif (type_1.parameters[i] <: AbstractDatum) && (has_ellipse(type_2) == Val{true})

            # can infer the datum    
            insert_param[i] = :(get_datum(X))
        end
    end            

    # and build the code
    if (length(type_1.parameters) == 1)
        quote
            add_param{X <: $(type_2)}(::Type{$(type_1)}, ::Type{X})            = $(type_1){$(insert_param[1])}    # it was not supplied
            add_param{T, X <: $(type_2)}(::Type{$(type_1){T}}, ::Type{X})      = $(type_1){T}                     # it was already supplied
            add_param{T <: $(type_1), X <: $(type_2)}(::Type{T}, ::X)          = add_param(T, X)
        end
    elseif (length(type_1.parameters) == 2)
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
# add basic maths to the point type
#####################################################

function add_maths{geodesy_type}(::Type{geodesy_type}, fields)

    nf = length(fields)  # number of fields and their names for the Geodesy variables

    # allow interaction with these types
    interact_types = [Vec{nf}, Vector]  

    # create an expression to add to
    q = quote end

    # and build (N.B. template functions keep clashing with FixedSizeArrays, so dont use them)
    for iT in interact_types

        if (nf == 2)
            qn = quote
                
                # addition 
                +(X::$(geodesy_type), dX::$(iT)) = typeof(X)(X.$(fields[1]) + dX[1], X.$(fields[2]) + dX[2])
                # +(dX::$(iT), X::$(geodesy_type)) = typeof(X)(X.$(fields[1]) + dX[1], X.$(fields[2]) + dX[2])  # should this version be a thing? Cant avoid a clash with FSA anyway...

                # subtraction
                -(X::$(geodesy_type), dX::$(iT)) = typeof(X)(X.$(fields[1]) - dX[1], X.$(fields[2]) - dX[2])

                # NaN
                isnan(X::$(geodesy_type)) = Vec{2, Bool}(isnan(X.$(fields[1])),   isnan(X.$(fields[2])))
                
            
            end
        elseif (nf == 3)
            qn = quote

                # addition 
                +(X::$(geodesy_type), dX::$(iT)) = typeof(X)(X.$(fields[1]) + dX[1], X.$(fields[2]) + dX[2], X.$(fields[3]) + dX[3])
                # +(dX::$(iT), X::$(geodesy_type)) = typeof(X)(X.$(fields[1]) + dX[1], X.$(fields[2]) + dX[2], X.$(fields[3]) + dX[3])   # should this version be a thing? Cant avoid a clash with FSA anyway...

                # subtraction
                -(X::$(geodesy_type), dX::$(iT)) = typeof(X)(X.$(fields[1]) - dX[1], X.$(fields[2]) - dX[2], X.$(fields[3]) - dX[3])

                # NaN
                isnan(X::$(geodesy_type)) = Vec{3, Bool}(isnan(X.$(fields[1])), isnan(X.$(fields[2])), isnan(X.$(fields[3])))
            end
        else
            qn = quote end
        end
        append!(q.args, qn.args)  # include them
    end
    
    return q

end


#####################################################
# Make the constructors force typing
#####################################################

function add_constructors{geodesy_type}(::Type{geodesy_type})

    q = quote

        # copy constructor
        call{T <: $(geodesy_type)}(::Type{T}, X::T) = X
        call{T <: $(geodesy_type)}(::Type{$(geodesy_type)}, X::T) = X
        call{T <: $(geodesy_type), U <: $(geodesy_type)}(::Type{T}, X::U) = geotransform(T, X, get_handler($(geodesy_type)), get_handler($(geodesy_type)))

    end

    # construction from scalars
    nfields = length(fieldnames(geodesy_type))
    if (nfields == 2)
        qn = quote

            # construction from scalars
            call(::Type{$(geodesy_type)}, x::Real, y::Real)                  = add_param($(geodesy_type))(x,y)
            call{T <: $(geodesy_type)}(::Type{T}, x::Real, y::Real, z::Real) = add_param(T)(x,y)   # legacy support, drop the height term
   
        end
    elseif (nfields == 3)
        
        qn = quote

            # construction from scalars
            call(::Type{$(geodesy_type)}, x::Real, y::Real, z::Real)         = add_param($(geodesy_type))(x,y,z)  
            call{T <: $(geodesy_type)}(::Type{T}, x::Real, y::Real)          = add_param(T)(x,y,0) # legacy support, insert the height term

        end
    else
        qn = quote
        #    call{T <: $(geodesy_type)}(::Type{T}, X...) = add_param(T)(X...)    # this can't be fast...
        end
    end 
    append!(q.args, qn.args)
    return q  
     
end 



#############################################################
# Add the cross constructors (constructing from other types)
#############################################################

function add_cross_constructors{type_1, type_2}(::Type{type_1}, ::Type{type_2})

    #h1, h2 = get_handler(type_1), get_handler(type_2)
    quote
        
        # no reference point
        call{T <: $(type_1)}(::Type{T}, X::$(type_2)) = geotransform(T, X)           # $(h1), $(h2)

        # with reference point
        call{T <: $(type_1)}(::Type{T}, X::$(type_2), ref) = geotransform(T, X, ref) # $(h1), $(h2)
    end
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
# construction
#
function add_vector_construction{geodesy_type}(::Type{geodesy_type}) 

    nfields = length(fieldnames(geodesy_type))
    if (nfields == 2)
        quote

            # allow construction from a matrix via convert as its value preserving
            function convert{T <: $(geodesy_type)}(::Type{Vector{T}}, X::Union{AbstractMatrix, Mat}; row::Bool=true)
                oT = add_param(T)
                n = (row) ? size(X,1) : size(X,2)
                Xout = Vector{oT}(n)  # cant make list comprehesion get the output type right
                if row 
                    for i = 1:n; Xout[i] = oT(X[i,1], X[i,2]); end
                else
                    for i = 1:n; Xout[i] = oT(X[1,i], X[2,i]); end
                end
                return Xout
            end
        end
    elseif (nfields == 3)
        quote
            # allow construction from a matrix via convert as its value preserving
            function convert{T <: $(geodesy_type)}(::Type{Vector{T}}, X::Union{AbstractMatrix, Mat}; row::Bool=true)
                oT = add_param(T)
                n = (row) ? size(X,1) : size(X,2)
                Xout = Vector{oT}(n)  # cant make list comprehesion get the output type right
                if (row)
                    for i = 1:n; Xout[i] = oT(X[i,1], X[i,2], X[i,3]); end
                else
                    for i = 1:n; Xout[i] = oT(X[1,i], X[2,i], X[3,i]); end
                end
                return Xout
            end
        end
    end
end




#
# Function build a code block with everything from above
# The idea is other packages can use this to create code blocks,
# and evaluate it in the Geodesy module to save exporting everything.  IDK
# 
function build_methods{geo_type}(::Type{geo_type}, fields=fieldnames(geo_type), known_types=[])

    # create a code block to amalgamate all the othe code blocks
    qb = Expr(:block)

    # add accessor functions
    append!(qb.args, Geodesy.add_accessors(geo_type, fields).args)
    
    # add the parameter manipulation methods
    append!(qb.args, Geodesy.add_param_methods(geo_type).args)

    # add basic maths
    append!(qb.args, Geodesy.add_maths(geo_type, fields).args)

    # add additional constructors
    append!(qb.args, Geodesy.add_constructors(geo_type).args)

    # the number of elements in the representation
    append!(qb.args, Geodesy.add_projection(geo_type).args)

    # add constructing vectors of these things from matrices
    append!(qb.args, Geodesy.add_vector_construction(geo_type).args)

    # NaN checks
    # append!(qb.args, Geodesy.add_nan_check(geo_type).args)

    # and add the binary methods
    for other_type in known_types

        # build both ways
        for comb in collect(permutations([geo_type, other_type]))
        
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








   
