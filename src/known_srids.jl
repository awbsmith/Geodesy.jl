
# list of srids for known datums (allows converting LLA / ECEF to SRID points)
# TODO: Check if these are actually correct!

# by default we don't know
SRID{T <: WorldPosition}(::T) = error("No known SRID / datum for a $(T)")
SRID{T <: WorldPosition}(::Type{T}) = error("No known SRID / datum for a $(T)")


############
# WGS84
############

# LLA
SRID(::Type{LLA_WGS84}) = SRID{:EPSG, 4326}        # EPSG code for lon lat wgs84 (GPS).  
SRID(::LLA_WGS84) = SRID(LLA_WGS84)

# ECEF
SRID(::Type{ECEF_WGS84}) = SRID{:EPSG, 4978}       # EPSG code for ecef wgs84 (GPS).  
SRID(::ECEF_WGS84) = SRID(ECEF_WGS84)


##############
# NAD27
##############

SRID(::Type{LLA{NAD27}}) = SRID{:EPSG, 4267}        
SRID(::LLA{NAD27}) = SRID(LLA{NAD27})

############## 
# ED50
##############

SRID(::Type{LLA{ED50}}) = SRID{:EPSG, 4230}        
SRID(::LLA{ED50}) = SRID(LLA{ED50})

############## 
# OSGB36
##############

# LLA
SRID(::Type{LLA{OSGB36}}) = SRID{:EPSG, 4277}        
SRID(::LLA{OSGB36}) = SRID(LLA{OSGB36})


############## 
# GDA94
##############

# LLA
SRID(::Type{LLA{GDA94}}) = SRID{:EPSG, 4283}        
SRID(::LLA{GDA94}) = SRID(LLA{GDA94})

############## 
# ETRS89
##############

# Europia
SRID(::Type{LLA{ETRS89}}) = SRID{:EPSG, 4258}        
SRID(::LLA{ETRS89}) = SRID(LA{ETRS89})


############## 
# NAD83
##############

SRID(::Type{LLA{NAD83}}) = SRID{:EPSG, 4269}        
SRID(::LLA{NAD83}) = SRID(LLA{NAD83})


