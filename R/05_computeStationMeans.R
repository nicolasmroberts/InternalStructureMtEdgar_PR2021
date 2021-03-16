# a list of station names for each core
stationList <- amsData$station

# Remove duplicate station names
stationList <- unique(stationList, incomparables = FALSE)

amsDataMeans = computeStationMeans(amsData, stationList)

# nameLength <- 9
# #this calls the locationPairing function and apppends location data to the amsDataset. 
# amsDataMeans <- locationPairing(amsDataMeans, locData, nameLength)
amsDataMeans$Name = amsDataMeans$station

# Convert data frame to simple feature object
# amsDataMeans$easting = as.numeric(amsDataMeans$easting) # make sure the coordinates are numerics (they default to character class)
# amsDataMeans$northing = as.numeric(amsDataMeans$northing)  # make sure the coordinates are numerics (they default to character class)
# amsDataMeans$elevation = as.numeric(replicate(length(amsDataMeans$easting), 10000))

# preload the two crs
crsDecDeg = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" 
crsUTM51S = "+proj=utm +zone=51 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

# convert amsDataMeans into an sf (simple feature) object
amsDataMeans = st_as_sf(amsDataMeans, coords = c("easting","northing"), crs = crsUTM51S)
