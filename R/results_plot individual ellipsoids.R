
# Plot ellipsoid for a single station
plotEllipsoid(amsDataMeans[amsDataMeans$station.x == "AME16_172",]$rotation, newAs[137], ellColor = "#dddddd")


ams_geom = st_geometry(amsDataMeans)

newLogAs = lapply(amsDataMeans$logA, function(s) normalizeOct(s,0.5))
newAs = lapply(newLogAs, function(s) exp(s))

plotEllipsoid(amsDataMeans$rotation, newAs, centers = ams_geom, numNonAdapt = 4, ellColor = "grey")



