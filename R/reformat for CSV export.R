# reformat amsDataMeans for exporting as a CSV file

amsDataMeans = st_transform(amsDataMeans, crs  = crsDecDeg)



xy <- do.call(rbind, st_geometry(amsDataMeans)) %>% 
     data.frame() %>% setNames(c("easting","northing"))

amsDataMeans$easting = xy$easting
amsDataMeans$northing = xy$northing

ellipsoids = lapply(1:length(amsDataMeans$Name.x), function(s) ellEllipsoidFromRotationLogA(amsDataMeans$rotation[[s]], amsDataMeans$logA[[s]]))

amsDataMeans$K11 = sapply(amsDataMeans$tensor, function(s) s[1,1])
amsDataMeans$K22 = sapply(amsDataMeans$tensor, function(s) s[2,2])
amsDataMeans$K33 = sapply(amsDataMeans$tensor, function(s) s[3,3])
amsDataMeans$K12 = sapply(amsDataMeans$tensor, function(s) s[1,2])
amsDataMeans$K13 = sapply(amsDataMeans$tensor, function(s) s[1,3])
amsDataMeans$K23 = sapply(amsDataMeans$tensor, function(s) s[2,3])

amsDataMeans$a1 = sapply(amsDataMeans$a, function(s) s[1])
amsDataMeans$a2 = sapply(amsDataMeans$a, function(s) s[2])
amsDataMeans$a3 = sapply(amsDataMeans$a, function(s) s[3])

k1 = data.frame(t(sapply(1:length(amsDataMeans$Name.x), function(i) geoTrendPlungeDegFromCartesian(lower(amsDataMeans$rotation[[i]][1,])))))
names(k1) = c("trend", "plunge")

k2 = data.frame(t(sapply(1:length(amsDataMeans$Name.x), function(i) geoTrendPlungeDegFromCartesian(lower(amsDataMeans$rotation[[i]][2,])))))
names(k2) = c("trend", "plunge")

k3 = data.frame(t(sapply(1:length(amsDataMeans$Name.x), function(i) geoTrendPlungeDegFromCartesian(lower(amsDataMeans$rotation[[i]][3,])))))
names(k3) = c("trend", "plunge")

amsDataMeans$K1dec = k1$trend
amsDataMeans$K1inc = k1$plunge
amsDataMeans$K2dec = k2$trend
amsDataMeans$K2inc = k2$plunge
amsDataMeans$K3dec = k3$trend
amsDataMeans$K3inc = k3$plunge

amsDataMeans$Pj = sapply(1:length(amsDataMeans$Name.x), function(i) ellJelinekP(amsDataMeans$logA[[i]]))
amsDataMeans$T = sapply(1:length(amsDataMeans$Name.x), function(i) ellLodeNu(amsDataMeans$logA[[i]]))

amsDataMeansForExport = amsDataMeans[,c("station.x",
                "easting",
                "northing",
                "specimens.x", 
                "Km","Pj","T",
                "K1dec", "K1inc",
                "K2dec", "K2inc",
                "K3dec", "K3inc",
                "K11","K22","K33","K12","K13","K23",
                "a1","a2","a3",
                "SUPERSUITE", "UNITNAME") ]

names(amsDataMeansForExport) = c("station",
                                 "easting",
                                 "northing",
                                 "specimens", 
                                 "Km","Pj","T",
                                 "K1dec", "K1inc",
                                 "K2dec", "K2inc",
                                 "K3dec", "K3inc",
                                 "K11","K22","K33","K12","K13","K23",
                                 "a1","a2","a3",
                                 "Supersuite", "Unit Name","geometry")

st_write(amsDataMeansForExport, "Appendices/amsDataMeansTable.csv")

# Now do the same thing for AMS specimen data

amsDataForExport = st_transform(amsData, crs  = crsDecDeg)

xy <- do.call(rbind, st_geometry(amsDataForExport)) %>% 
     data.frame() %>% setNames(c("easting","northing"))

amsDataForExport$easting = xy$easting
amsDataForExport$northing = xy$northing

amsDataForExport = amsDataForExport[,c("Name",
                                        "easting",
                                        "northing",
                                        "Km","Pj","T",
                                        "K1dec", "K1inc",
                                        "K2dec", "K2inc",
                                        "K3dec", "K3inc",
                                        "K11","K22","K33","K12","K13","K23",
                                        "a1","a2","a3",
                                        "SUPERSUITE", "UNITNAME") ]

names(amsDataForExport) = c("core name",
                            "easting",
                            "northing",
                            "Km","Pj","T",
                            "K1dec", "K1inc",
                            "K2dec", "K2inc",
                            "K3dec", "K3inc",
                            "K11","K22","K33","K12","K13","K23",
                            "a1","a2","a3",
                            "Supersuite", "Unit Name","geometry")

st_write(amsDataForExport, "Appendices/amsDataTable.csv")

# Hysteresis and Thermalsusc table

magMinTableForExport = st_transform(magMinLocations, crs  = crsDecDeg)

xy <- do.call(rbind, st_geometry(magMinTableForExport)) %>% 
     data.frame() %>% setNames(c("easting","northing"))

magMinTableForExport$easting = xy$easting
magMinTableForExport$northing = xy$northing

magMinTableForExport = magMinTableForExport[,c("Name",
                                       "easting",
                                       "northing", "VSM","ThermSusc",
                                       "SUPERSUITE", "UNITNAME") ]

names(magMinTableForExport) = c("Name",
                                "easting",
                                "northing", "VSM","ThermSusc",
                                "SUPERSUITE", "UNITNAME","geometry")

st_write(magMinTableForExport, "Appendices/amsMagMinTable.csv")



#FieldFabrics

# Hysteresis and Thermalsusc table

fieldFabGranites2016ForExport = st_transform(fieldFabGranites2016, crs  = crsDecDeg)
fieldFabGranites2016ForExport = fieldFabGranites2016ForExport[,c("Name","strike","dip","Rake1","SUPERSUITE","UNITNAME")]

names(fieldFabGranites2016ForExport) = c("station","strike","dip","rake","Supersuite","Unit Name","geometry")


fieldFabGranites2018ForExport = st_transform(fieldFabGranites2018, crs  = crsDecDeg)
fieldFabGranites2018ForExport = fieldFabGranites2018ForExport[,c("fieldStation","strike","dip","assocRake","SUPERSUITE","UNITNAME")]
names(fieldFabGranites2018ForExport) = c("station","strike","dip","rake","Supersuite","Unit Name","geometry")



fieldFabForExport = rbind(fieldFabGranites2016ForExport,fieldFabGranites2018ForExport)


xy <- do.call(rbind, st_geometry(fieldFabForExport)) %>% 
     data.frame() %>% setNames(c("easting","northing"))

fieldFabForExport$easting = xy$easting
fieldFabForExport$northing = xy$northing

fieldFabForExport = fieldFabForExport[,c("station","easting","northing", "strike","dip","rake","Supersuite","Unit Name") ]


st_write(fieldFabForExport, "Appendices/fieldFabrics.csv")






