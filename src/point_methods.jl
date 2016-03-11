
###########################################
# defines methods to add to point types
# Note the "geotransform" function handles all transformations
# while the "convert" function converts types leaving values unchanged
###########################################

# TODO: This whole thing should probably be switched to a codegen format

# Convenient union types 
Geodesy_fam = Union{WorldPosition, WorldSurfacePosition, LocalPosition}  # try to ctach everything
World_fam = Union{WorldPosition, WorldSurfacePosition}
Local_fam = Union{LocalPosition}
LL_fam = Union{LLA, LL}


Proj4_fam = Union{WorldPosition, Proj4Type}      # acceptable types to give to Proj4

Vec3_fam = Union{WorldPosition, LocalPosition}   # for three element point types
Vec2_fam = Union{LL}                             # for three element point types

ELL_param_fam = Union{LLA, LL, ECEF}             # things where an ellipse / psuedo datum are the template param
LL_param_fam = Union{ENU}                        # things where an LLA points is the template param


#########################################
# Accessors
#########################################

# Point translators (lgacy support)
getX{T <: Geodesy_fam}(X::T) = X[1]
getY{T <: Geodesy_fam}(X::T) = X[2]
getZ{T <: Vec3_fam}(X::T) = X[3]

get_lon{T <: Geodesy_fam}(X::T) = X[1]
get_lat{T <: Geodesy_fam}(X::T) = X[2]
get_alt{T <: Vec3_fam}(X::T) = X[3]

get_east{T <: Geodesy_fam}(X::T) = X[1]
get_north{T <: Geodesy_fam}(X::T) = X[2]
get_up{T <: Vec3_fam}(X::T) = X[3]





#####################################################
# Template parameter manipultation
#####################################################

# replace a TypeVar parameter with a DataType parameter where needed

#TODO: Should these be promote rules or something?

#add_param{T <: ELL_param_fam}(::Type{T}) = typeof(T.parameters[1]) == DataType ? T : T{UnknownDatum}  # not type safe :-(
@generated add_param{T <: ELL_param_fam}(::Type{T}) = (typeof(T.parameters[1]) == DataType) ? :(T) :  :(T{UnknownDatum})

#add_param{T <: LL_param_fam}(::Type{T}) = typeof(T.parameters[1]) == DataType ? T : T{UnknownRef}      # not type safe :-(
@generated add_param{T <: LL_param_fam}(::Type{T}) = (typeof(T.parameters[1]) == DataType) ? :(T) :  :(T{UnknownRef})

add_param{T <: CRS}(::Type{T}) = typeof(T.parameters[1]) == DataType ? T : error("always specify an SRID when using the CRS type")

# retrieve datums and ellipsoids
@generated ellipsoid{T <: ELL_param_fam}(::Type{T}) = :($(ellipsoid(T.parameters[1])))           # reference ellipsoid for the position

# methods to pull the reference point from the template
@generated ELL_type{T <: ELL_param_fam}(::Type{T}) = :($(T.parameters[1]))
@generated ELL_type{T <: ELL_param_fam}(::T) = :($(T.parameters[1]))
@generated ELL_type{T <: LL_param_fam}(::Type{T}) = :($(ELL_type(T.parameters[1])))
@generated LL_ref{T <: LL_param_fam}(::Type{T}) = :($(T.parameters[1]))
@generated LL_ref{T <: LL_param_fam}(::T) = :($(T.parameters[1]))

# add the ENU ref or chnage it to LL if its LLA
@generated add_LL_ref{T <: ENU}(::Type{T}) = typeof(T.parameters[1]) == TypeVar ? :(ENU_NULL) : ((typeof(T.parameters[1]) <: LLA) ? :($(ENU{LL(T.parameters[1])})) :  :($(ENU{T.parameters[1]})))


#####################################################
# Force the template parameter when constructing
#####################################################

call(::Type{LL}, x::Real, y::Real) = add_param(LL)(x,y)
call(::Type{LLA}, x::Real, y::Real, z::Real) = add_param(LLA)(x,y,z)
call(::Type{ECEF}, x::Real, y::Real, z::Real) = add_param(ECEF)(x,y,z)
call(::Type{ENU}, x::Real, y::Real, z::Real) = add_param(ENU)(x,y,z)


#
# allow crazy parameter construction for ALL types (legacy support)
# 
call{T <: Vec3_fam}(::Type{T}, x::Real, y::Real) = add_param(T)(x,y,0.0)
call{T <: Vec2_fam}(::Type{T}, x::Real, y::Real, z::Real) = add_param(T)(x,y)



####################################
### Proj4 projections
### Try to get the SRID for the point type
####################################

# build methods to get a proj 4 projection (a coordinate reference system) for stuff we know about
get_projection{T <: Union{WorldPosition, WorldSurfacePosition}}(X::T) = get_projection(T)
get_projection{T <: Union{WorldPosition, WorldSurfacePosition}}(::Type{T}) = get_projection(SRID(T))



#####################################################
# Build abstract calls for the transformations
#####################################################


#
# Allow LL <-> LLA
#
geotransform{T}(::Type{LL{T}}, lla::LLA{T}) = LL{T}(lla.lon, lla.lat)
geotransform{T}(::Type{LL}, lla::LLA{T}) = LL{T}(lla.lon, lla.lat)

geotransform{T}(::Type{LLA{T}}, ll::LL{T}) = LLA{T}(lla.lon, lla.lat, 0.0)
geotransform{T}(::Type{LLA}, ll::LL{T}) = LLA{T}(lla.lon, lla.lat, 0.0)


# No conversions, stripping or adding the template parameters
macro default_convs(type_name, known_subtype, null_subtype)
    eval(quote

        # geotransform it to itself (i.e. do nothing)
        geotransform{T}(::Type{$(type_name){T}}, X::$(type_name){T}) = X

        # allow stripping the template parameter from the template
        geotransform{T <: $(known_subtype)}(::Type{$(type_name){$(null_subtype)}}, X::$(type_name){T}) = $(type_name){$(null_subtype)}(X...)

        # allow inserting a template parameter into the template
        geotransform{T <: $(known_subtype)}(::Type{$(type_name){T}}, X::$(type_name){$(null_subtype)}) = $(type_name){T}(X...)

    end)
end

# generate default conversions for everything
@default_convs(LLA, KnownDatum, UnknownDatum)
@default_convs(LL, KnownDatum, UnknownDatum)
@default_convs(ECEF, KnownDatum, UnknownDatum)
@default_convs(ENU, LLA, UnknownRef)

# dont macro srid type's, we don't want to allow stripping of the SRID 
geotransform{T <: SRID}(::Type{CRS}, X::CRS{T}) = X  
geotransform{T <: SRID}(::Type{CRS{T}}, X::CRS{T}) = X               


#
# templates for constructing one point type from another (I want to leave convert free to do value preserving conversions)
#
macro add_cross_const(Type1, Type2)
    eval(quote

        # copy constructor for the first one
        call{T <: $(Type1), U <: $(Type1)}(::Type{T}, X::U) = geotransform(T, X)

        # and the cross versions
        call{T <: $(Type1), U <: $(Type2)}(::Type{T}, X::U) = geotransform(T, X)
        call{T <: $(Type2), U <: $(Type1)}(::Type{T}, X::U) = geotransform(T, X)
    end)
end
@add_cross_const(WorldPosition, LocalPosition)
@add_cross_const(LocalPosition, WorldSurfacePosition)
@add_cross_const(WorldSurfacePosition, WorldPosition)


#
# constructing one point type from another and a reference
#
call{T <: Local_fam}(::Type{T}, X::World_fam, ll_ref::LL_fam) = geotransform(T, X, ll_ref)
call{T <: World_fam}(::Type{T}, X::Local_fam, ll_ref::LL_fam) = geotransform(T, X, ll_ref)


#
# allow construction from a matrix via convert as its value preserving
#
function convert{T <: Geodesy_fam}(::Type{Vector{T}}, X::AbstractMatrix; row::Bool=true)
    oT = add_param(T)
    n = (row) ? size(X,1) : size(X,2)
    Xout = Vector{oT}(n)  # cant make list comprehesion get the output type right
    if (T <: Vec2_fam) && row 
        for i = 1:n; Xout[i] = oT(X[i,1], X[i,2]); end
    elseif (row)
        for i = 1:n; Xout[i] = oT(X[i,1], X[i,2], X[i,3]); end
    elseif (T <: Vec2_fam)
        for i = 1:n; Xout[i] = oT(X[1,i], X[2,i]); end
    else
        for i = 1:n; Xout[i] = oT(X[1,i], X[2,i], X[3,i]); end
    end
    return Xout
end





###################################################################
### Methods to add template parameters based on other arguments ###
###################################################################

# we want to make sure any created Vector have the template parameter in them
# add_param{T <: ELL_param_fam, U <: ELL_param_fam}(::Type{T}, ::Type{U}) = T ==  add_param(T) ? T : T{U.parameters[1]}                              # not type safe :-(
@generated add_param{T <: ELL_param_fam, U <: ELL_param_fam}(::Type{T}, ::Type{U}) = (T == add_param(T)) ? :(T) : :(T{$(ELL_type(U))})

# add_param{T <: ELL_param_fam, U <: LL_param_fam}(::Type{T}, ::Type{U}) = T ==  add_param(T) ? T : T{typeof(U.parameters[1]).parameters[1]}        # not type safe :-(
@generated add_param{T <: ELL_param_fam, U <: LL_param_fam}(::Type{T}, ::Type{U}) = T ==  add_param(T) ? :(T) : :(T{$(typeof(U.parameters[1]).parameters[1])})

# default to not include reference position in local types
add_param{T <: LL_param_fam, U <: Geodesy_fam}(::Type{T}, ::Type{U}) = add_param(T)  

# add parameters when it might be an SRID
@generated add_param{T <: Geodesy_fam, U <: Geodesy_fam}(::Type{T}, ::Type{U}) = (T <: CRS) ? :(T) : ((U <: CRS) ? :(add_param(T)) : :(add_param(T, U)))




####################################
### Transformations of Vectors   ###
####################################

# a vectorized way to perform transformations
function geotransform{T <: Geodesy_fam, U <: Geodesy_fam}(::Type{T}, X::Vector{U})

    # make sure the output parameter is filled
    oT =  add_param(T, U)

    # is Proj4 involved?
    if (T <: CRS) || (U <: CRS)
        Xout = proj4_vectorized(oT, X)
    else
        Xout = Vector{oT}(length(X))
        @inbounds for i = 1:length(X)
            Xout[i] = geotransform(oT, X[i])
        end
    end
    return Xout
end


# a vectorized way to perform transformations with reference points
function geotransform{T <: Geodesy_fam, U <: Geodesy_fam, V <: LL_fam}(::Type{T}, X::Vector{U}, ll_ref::V)

    # get inputs and desired output types
    oT = add_param(T, V)

    Xout = Vector{oT}(length(X))
    @inbounds for i = 1:length(X)
        Xout[i] = geotransform(oT, X[i], ll_ref)
    end
    return Xout
end


# worker code for proj4 vectorization (looping prj4 is slow)
function proj4_vectorized{T <: Proj4_fam, U <: Proj4_fam}(::Type{T}, X::Vector{U})

    if !((T <: CRS) || (U <: CRS))
        warn("Unexpected: using Proj4 to geotransform between Geodesy point types.  How'd this happen")
    end

    # convert to a matrix 
    mat = Matrix{Float64}(length(X), 3)
    @inbounds for i = 1:length(X)
        mat[i, 1] = X[i][1]
        mat[i, 2] = X[i][2]
        mat[i, 3] = U <: Vec2_fam ? 0.0 : X[i][3]
    end

    # perform it
    Proj4.transform!(Geodesy.get_projection(U), Geodesy.get_projection(T), mat, false)

    # and assign the output
    X = Vector{T}(length(X))
    @inbounds for i = 1:length(X)
        X[i] = T(mat[i,1], mat[i,2], mat[i,3])
    end
    
    return X
end







   
