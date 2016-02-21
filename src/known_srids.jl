
# list of srids for known datums (allows converting LLA / ECEF to SRID points)
# TODO: Check if these are actually correct!

# by default we don't know
SRID{T <: WorldPosition}(::Type{T}) = error("No known SRID / datum for a $(T)")
SRID{T <: WorldPosition}(::T) = SRID(T)


############
# WGS84
############

# LLA
#SRID(::Type{LLA_WGS84}) = SRID{:EPSG, 4326}        # EPSG code for lon lat wgs84 (GPS).  

# ECEF
SRID(::Type{ECEF_WGS84}) = SRID{:EPSG, 4978}       # EPSG code for ecef wgs84 (GPS).  


##############
# NAD27
##############

SRID(::Type{LLA{NAD27}}) = SRID{:EPSG, 4267}        

############## 
# ED50
##############

SRID(::Type{LLA{ED50}}) = SRID{:EPSG, 4230}        

############## 
# OSGB36
##############

# LLA
SRID(::Type{LLA{OSGB36}}) = SRID{:EPSG, 4277}        


############## 
# GDA94
##############

# LLA
SRID(::Type{LLA{GDA94}}) = SRID{:EPSG, 4283}        

############## 
# ETRS89
##############

# Europia
SRID(::Type{LLA{ETRS89}}) = SRID{:EPSG, 4258}        


############## 
# NAD83
##############

SRID(::Type{LLA{NAD83}}) = SRID{:EPSG, 4269}        


