
###########################################
# defines methods to add to point types
# Note the "transform" function handles all transformations
# while the "convert" function converts types leaving values unchanged
###########################################

# Convenient union types 
Geodesy_fam = Union{WorldPosition, WorldSurfacePosition, LocalPosition}  # try to ctach everything
World_fam = Union{WorldPosition, WorldSurfacePosition}
Local_fam = Union{LocalPosition}
LL_fam = Union{LLA, LL}

Proj4_fam = Union{WorldPosition, SRID_Pos}          # acceptable types to give to Proj4

Vec3_fam = Union{WorldPosition, LocalPosition}  # for three element point types
Vec2_fam = Union{LL}  # for three element point types

ELL_param_fam = Union{LLA, LL, ECEF}  # things where an ellipse / psuedo datum are the template param
LLA_param_fam = Union{ENU}            # things where an LLA points is the template param


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

# replace TypeVar parameter with a DataType parameter where needed
add_param{T <: ELL_param_fam}(::Type{T}) = typeof(T.parameters[1]) == DataType ? T : T{UnknownEllipse}
add_param{T <: LLA_param_fam}(::Type{T}) = typeof(T.parameters[1]) == DataType ? T : T{UnknownRef}
add_param{T <: SRID_Pos}(::Type{T}) = typeof(T.parameters[1]) == DataType ? T : error("always specify an SRID when using the SRID_Pos type")

# retrieve datums and ellipsoids
ellipsoid{T <: ELL_param_fam}(::Type{T}) = ellipsoid(T.parameters[1])           # reference ellipsoid for the position

# methods to pull the reference point from the template
LLA_ref{T <: LLA_param_fam}(::Type{T}) = T.parameters[1]
LLA_ref{T <: LLA_param_fam}(::T) = T.parameters[1]
ELL_type{T <: ELL_param_fam}(::Type{T}) = T.parameters[1]
ELL_type{T <: ELL_param_fam}(::T) = T.parameters[1]


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
get_projection{T <: Union{WorldPosition, WorldSurfacePosition}}(X::T) = get_projection(SRID(T))
get_projection{T <: Union{WorldPosition, WorldSurfacePosition}}(::Type{T}) = get_projection(SRID(T))



#####################################################
# Build abstract calls for the transformations
#####################################################


#
# Allow LL <-> LLA
#
transform{T}(::Type{LL{T}}, lla::LLA{T}) = LL{T}(lla.lon, lla.lat)
transform{T}(::Type{LL}, lla::LLA{T}) = LL{T}(lla.lon, lla.lat)

transform{T}(::Type{LLA{T}}, ll::LL{T}) = LLA{T}(lla.lon, lla.lat, 0.0)
transform{T}(::Type{LLA}, ll::LL{T}) = LLA{T}(lla.lon, lla.lat, 0.0)


# No conversions, stripping or adding the template parameters
macro default_convs(type_name, known_subtype, null_subtype)
	eval(quote

		# transform it to itself (i.e. do nothing)
		transform{T}(::Type{$(type_name){T}}, X::$(type_name){T}) = X

		# allow stripping the template parameter from the template
		transform{T <: $(known_subtype)}(::Type{$(type_name){$(null_subtype)}}, X::$(type_name){T}) = $(type_name){$(null_subtype)}(X...)

		# allow inserting a template parameter into the template
		transform{T <: $(known_subtype)}(::Type{$(type_name){T}}, X::$(type_name){$(null_subtype)}) = $(type_name){T}(X...)

    end)
end

# generate default conversions for everything
@default_convs(LLA, KnownEllipse, UnknownEllipse)
@default_convs(LL, KnownEllipse, UnknownEllipse)
@default_convs(ECEF, KnownEllipse, UnknownEllipse)
@default_convs(ENU, LLA, UnknownRef)

# dont macro srid type's, we don't want to allow stripping of the SRID 
transform{T <: SRID}(::Type{SRID_Pos}, X::SRID_Pos{T}) = X  
transform{T <: SRID}(::Type{SRID_Pos{T}}, X::SRID_Pos{T}) = X		   	


#
# templates for constructing one point type from another (I want to leave convert free to do value preserving conversions)
#
macro add_cross_const(Type1, Type2)
	eval(quote

		# copy constructor for the first one
		call{T <: $(Type1), U <: $(Type1)}(::Type{T}, X::U) = transform(T, X)

		# and the cross versions
		call{T <: $(Type1), U <: $(Type2)}(::Type{T}, X::U) = transform(T, X)
		call{T <: $(Type2), U <: $(Type1)}(::Type{T}, X::U) = transform(T, X)
    end)
end
@add_cross_const(WorldPosition, LocalPosition)
@add_cross_const(LocalPosition, WorldSurfacePosition)
@add_cross_const(WorldSurfacePosition, WorldPosition)


#
# constructing one point type from another and a reference
#
call{T <: Local_fam}(::Type{T}, X::World_fam, lla_ref::LL_fam) = transform(T, X, lla_ref)
call{T <: World_fam}(::Type{T}, X::Local_fam, lla_ref::LL_fam) = transform(T, X, lla_ref)


#
# allow construction from a matrix
#
function call{T <: Geodesy_fam, U <: Union{Val{:col}, Val{:row}}}(::Type{T}, X::AbstractMatrix, ::Type{U})
	oT = add_param(T)
	n = (U == Val{:col}) ? size(X,2) : size(X,1)
	Xout = Vector{oT}(n)
	for i = 1:n
		Xout[i] = (U == Val{:col}) ? oT(X[1,i], X[2,i], X[3,i]) : oT(X[i,1], X[i,2], X[i,3])
	end
	return Xout
end





####################################
### Transforms on Vectors        ###
####################################

# we want to make sure any created Vector have the template parameter in them
add_param{T <: ELL_param_fam, U <: ELL_param_fam}(::Type{T}, ::Type{U}) = T ==  add_param(T) ? T : T{U.parameters[1]}
add_param{T <: ELL_param_fam, U <: LLA_param_fam}(::Type{T}, ::Type{U}) = T ==  add_param(T) ? T : T{typeof(U.parameters[1]).parameters[1]}
add_param{T <: LLA_param_fam, U <: Geodesy_fam}(::Type{T}, ::Type{U}) = add_param(T)  # default to not include reference position in local types


# a vectorized way to perform transformations
function transform{T <: Geodesy_fam, U <: Geodesy_fam}(::Type{T}, X::Vector{U})

	# make sure the output parameter is filled
	oT = (T <: SRID_Pos) ? T : ((U <: SRID_Pos) ? add_param(T) : oT = add_param(T, U))

	# is Proj4 involved?
	if (T <: SRID_Pos) || (U <: SRID_Pos)
		Xout = proj4_vectorized(oT, X)
	else
		Xout = Vector{oT}(length(X))
		@inbounds for (i, Xp) in enumerate(X)
			Xout[i] = transform(oT, Xp)
		end
	end
	return Xout
end


# a vectorized way to perform transformations with reference points
function transform{T <: Geodesy_fam, U <: Geodesy_fam, V <: LL_fam}(::Type{T}, X::Vector{U}, lla_ref::V)

	# get inputs and desired output types
	oT = add_param(T, V)

	Xout = Vector{oT}(length(X))
	@inbounds for (i, Xp) in enumerate(X)
		Xout[i] = transform(oT, Xp, lla_ref)
	end
	return Xout
end


# when converting vectors of srid points, do it with a single call to Proj4.transform
# I haven't been able to make repeated calls to Proj4.transform work
function proj4_vectorized{T <: Proj4_fam, U <: Proj4_fam}(::Type{T}, X::Vector{U})

	if !((T <: SRID_Pos) || (U <: SRID_Pos))
		warn("Unexpected: using Proj4 to transform between Geodesy point types.  How'd this happen")
	end

	# convert to a matrix 
	mat = Matrix{Float64}(length(X), 3)
	@inbounds for i = 1:length(X)
		mat[i, 1] = X[i][1]
		mat[i, 2] = X[i][2]
		mat[i, 3] = X[i][3]
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








   
