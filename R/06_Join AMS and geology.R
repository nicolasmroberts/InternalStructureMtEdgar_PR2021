# Perform a join of the geology shapefile data and the amsData--basically, add the unit, supersuite, etc. name to all AMS data points
amsDataMeans = st_join(amsDataMeans, mtEdgarMap, join = st_within)
amsDataMeans = st_join(amsDataMeans, ellipsoidTriaxiality, join = st_equals)

# Exclude data from the Corunna Downs dome
amsDataMeans = amsDataMeans[amsDataMeans$UNITNAME != "Carbana Monzogranite" & amsDataMeans$UNITNAME != "Apex Basalt",]

amsData = st_join(amsData, mtEdgarMap, join = st_within)

# Exclude data from the Corunna Downs dome
amsData = amsData[amsData$UNITNAME != "Carbana Monzogranite" & amsData$UNITNAME != "Apex Basalt" ,]

# Join geochron and geology
geochron = st_join(geochron, mtEdgarMap, join = st_within)

# Join other datasets
DayPlotData = st_join(DayPlotData, mtEdgarMap, join = st_within)
fieldFab2018 = st_join(fieldFab2018, mtEdgarMap, join = st_within)
fieldFabGranites2018 = fieldFab2018[is.na(fieldFab2018$SUPERSUITE) == FALSE,]
fieldFabGranites2018 =  fieldFabGranites2018[fieldFabGranites2018$UNITNAME != "Carbana Monzogranite" & fieldFabGranites2018$UNITNAME != "Apex Basalt",]

fieldFab2016 = st_join(fieldFab2016, mtEdgarMap, join = st_within)
fieldFabGranites2016 = fieldFab2016[is.na(fieldFab2016$SUPERSUITE) == FALSE,]
fieldFabCores = fieldFabGranites2016[fieldFabGranites2016$Core == 1,] 

microStructures = st_join(microStructures, mtEdgarMap, join = st_within)
microStructures = microStructures[microStructures$UNITNAME != "Carbana Monzogranite" & microStructures$UNITNAME != "Apex Basalt",]

magMinLocations = st_join(magMinLocations, mtEdgarMap, join = st_within)
magMinLocations = magMinLocations[magMinLocations$UNITNAME != "Carbana Monzogranite" & magMinLocations$UNITNAME != "Apex Basalt",]



                  