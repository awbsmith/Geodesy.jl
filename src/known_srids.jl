
# list of srids for known datums (allows converting LLA / ECEF to SRID points)
# TODO: Check if these are actually correct!

# by default we don't know
SRID{T}(::Type{T}) = error("No known SRID for type $(T).\nPlease overload the SRID constructor to return an SRID for the $(T) point type:\nSRID(::Type{$(T)}) = SRID{..., ...}")
SRID{T}(::T) = SRID(T)

# build the list
known_srids = [

############
# WGS84
############

LLA{WGS84}      :EPSG       4326;    # EPSG code for lon lat wgs84 (GPS).
ECEF{WGS84}     :EPSG       4978;    # EPSG code for lon ecef wgs84 (GPS).  

##############
# NAD27
##############

LLA{NAD27}      :EPSG       4267;


############## 
# ED50
##############

LLA{ED50}       :EPSG       4230;

############## 
# OSGB36
##############

LLA{OSGB36}     :EPSG       4277;

############## 
# GDA94
##############

# Australia
LLA{GDA94}     :EPSG        4283;


############## 
# ETRS89
##############

# Europia
LLA{ETRS89}    :EPSG        4258;


############## 
# NAD83
##############

# Americania
LLA{NAD83}     :EPSG        4269

]

# and code gen them up
for i = 1:size(known_srids,1)

    q = quote
    
        # constructor from type
        #SRID(::Type{$(known_srids[i,1])}) = SRID{$(known_srids[i,2]), $(known_srids[i,3])}
        SRID(::Type{$(known_srids[i,1])}) = SRID{:EPSG, $(known_srids[i,3])}         # want the above line but I'm bad at Julia

        # constructor from a point
        SRID(::$(known_srids[i,1])) = SRID($(known_srids[i,1]))

        # add a get_srid method (this should add it to the TrSRIDType trait).  Is this a good thing?)
        get_srid(::Type{$(known_srids[i,1])}) = SRID($(known_srids[i,1]))
        get_srid(::$(known_srids[i,1])) = SRID($(known_srids[i,1]))

    end
    eval(q)
end


        


