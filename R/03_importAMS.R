#Path to AGICO txt data file
AGICOcsvPath <- "data/amsData_12Sept_2019.csv"

#load the AGICO csv file
dataFrame <- geoDataFromFile(AGICOcsvPath, separator = ",")

#This computes useful measures from the raw AGICO csv. This is the main variable which all the statistics will draw from. 
amsData <- geoEllipsoidDataFromAGICOFile(AGICOcsvPath, separator = ",", sapply(dataFrame$Km, function(s) s*10^(-6)), doNormalize=TRUE)
#for some reason, AGICO stores all the bulk susceptibilities as the actual susceptibility * 10^6, so this quickly corrects it.
amsData$Km <- sapply(amsData$Km, function(s) s*10^(-6))

#the file path to the station name/location csv (see below for details on format)
locationPath <- "data/sampleLocations.csv"

#Loads a csv file with the location (easting, northing) information. The table needs to be formatted with headings "name", "easting", and "northing"--all lower case, no spaces. The table can include more stations than the amsData have, but must inlcude all stations in amsData. More details below. 
locData <- read.table(locationPath, sep = ",", header=TRUE,stringsAsFactors = FALSE)

nameLength <- 9

#this calls the locationPairing function and apppends location data to the amsDataset. 
amsData <- locationPairing(amsData, locData, nameLength)


# Convert data frame to simple feature object
amsData$easting = as.numeric(amsData$easting) # make sure the coordinates are numerics (they default to character class)
amsData$northing = as.numeric(amsData$northing)  # make sure the coordinates are numerics (they default to character class)
#amsData$elevation = as.numeric(replicate(length(amsData$easting), 10000))

# preload the two crs
crsDecDeg = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" 
crsUTM51S = "+proj=utm +zone=51 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

# convert amsData into an sf (simple feature) object
amsData = st_as_sf(amsData, coords = c("easting","northing"), crs = crsDecDeg)

# convert amsData from decimal degrees to UTM 51S meters
amsData = st_transform(amsData, crs  = crsUTM51S)

#----------------------
#
micro = read.csv("data/graniteMicrostructures.csv")
microStructures = micro[micro$Classification != "",]

microStructures$Classification = factor(microStructures$Classification, levels = c( "M", "SM2", "SM3", "HT"))

crsDecDeg = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" 
crsUTM51S = "+proj=utm +zone=51 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

# convert amsData into an sf (simple feature) object
microStructures = st_as_sf(microStructures, coords = c("easting","northing"), crs = crsDecDeg)

# convert amsData from decimal degrees to UTM 51S meters
microStructures = st_transform(microStructures, crs  = crsUTM51S)

#-----------------------------
magMinLocations = read.csv("data/Magnetics/MagMinNames.csv")

magMinLocations = locationPairing(magMinLocations, locData, 9)

crsDecDeg = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" 
crsUTM51S = "+proj=utm +zone=51 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

# convert amsData into an sf (simple feature) object
magMinLocations = st_as_sf(magMinLocations, coords = c("easting","northing"), crs = crsDecDeg)

# convert amsData from decimal degrees to UTM 51S meters
magMinLocations = st_transform(magMinLocations, crs  = crsUTM51S)

#-----------------------

ellipsoidTriaxiality = read.csv("data/stationCategorization_triaxiality.csv")

ellipsoidTriaxiality = locationPairing(ellipsoidTriaxiality, locData, 9)

crsDecDeg = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" 
crsUTM51S = "+proj=utm +zone=51 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

# convert amsData into an sf (simple feature) object
ellipsoidTriaxiality = st_as_sf(ellipsoidTriaxiality, coords = c("easting","northing"), crs = crsDecDeg)

# convert amsData from decimal degrees to UTM 51S meters
ellipsoidTriaxiality = st_transform(ellipsoidTriaxiality, crs  = crsUTM51S)


#------------------------
#
fieldFab2018 = read.csv("data/FieldData2018.csv")

fieldFab2018$pole = lapply(1:nrow(fieldFab2018), function(i) geoCartesianFromStrikeDipDeg(c(fieldFab2018$strike[i], fieldFab2018$dip[i] )))
# convert amsData into an sf (simple feature) object
# 
fieldFab2018 = st_as_sf(fieldFab2018, coords = c("easting","northing"), crs = crsDecDeg)

# convert amsData from decimal degrees to UTM 51S meters
fieldFab2018 = st_transform(fieldFab2018, crs  = crsUTM51S)

# 
fieldFab2016 = read.csv("data/FieldData2016.csv")

fieldFab2016$pole = lapply(1:nrow(fieldFab2016), function(i) geoCartesianFromStrikeDipDeg(c(fieldFab2016$strike[i], fieldFab2016$dip[i] )))
# convert amsData into an sf (simple feature) object
# 
fieldFab2016 = st_as_sf(fieldFab2016, coords = c("easting","northing"), crs = crsDecDeg)

# convert amsData from decimal degrees to UTM 51S meters
fieldFab2016 = st_transform(fieldFab2016, crs  = crsUTM51S)

