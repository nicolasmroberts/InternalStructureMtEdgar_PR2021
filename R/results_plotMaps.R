
# Plot the gamma ray spectrometry imagery
pctKMap = ggplot() +
     geom_raster(data = gammaKFrame, aes(x = x, y = y, fill = (UnfilteredPctK)), show.legend = FALSE) +
     scale_fill_gradient(low = "black", high = "white" , limits = c(0.5, 4.5), oob = squish) +
     coord_fixed(expand = FALSE) +
     xlim(UTMxlim) +
     ylim(UTMylim) +
     theme_bw() +
     theme(axis.title = element_blank())

# Plot geology with full color
geoMapMtEdgar = ggplot() +
     geom_sf(data = mtEdgarMap["UNITNAME"], aes(fill = UNITNAME), size = 0.25,  show.legend = FALSE) +
     scale_fill_manual(values = pilbaraColors) +
     theme_bw() +
     theme(axis.title = element_blank()) +
     coord_sf(expand = FALSE, datum = crsDecDeg)

# Plot granite geology in greyscale 
geoMapMtEdgarSuperSuitesBW = ggplot() +
  geom_sf(data = mtEdgarMap["SUPERSUITE"], aes(fill = SUPERSUITE), size = 0.25,  show.legend = FALSE, inherit.aes = FALSE) +
  scale_fill_manual(values = superSuiteColorsBW, guide = FALSE) +
  theme_bw() +
  theme(axis.title = element_blank()) +
  coord_sf(expand = FALSE, datum = crsDecDeg)


# plot map with geochron and ams locations
geochronMtEdgar = geochron[geochron$ANALYSIS_A > 3100,]
geochronMtEdgar$ANALYSIS_A = round(geochronMtEdgar$ANALYSIS_A/1000, digits = 2)
geochronMtEdgar = geochronMtEdgar[is.na(geochronMtEdgar$SUPERSUITE) == FALSE,]

xyGeochron <- do.call(rbind, st_geometry(geochronMtEdgar)) %>% 
  as_tibble() %>% setNames(c("x","y"))

geochronOnGeoMap = geoMapMtEdgar + 
  geom_sf(data = geochronMtEdgar, aes(),color = "black", fill = "white", shape = 21, size = 2) +
  annotate("rect", xmin = mtEdgarBboxUTM[[2]] - 21000, xmax = mtEdgarBboxUTM[[2]] - 1000, ymin = mtEdgarBboxUTM[[3]] + 1000, ymax = mtEdgarBboxUTM[[3]] + 4000, fill = "black") + 
  mapply(function(xx,yy, dates){annotate("text", x = xx + 1500, y = yy +1500, label = paste(dates), size = 3)},xyGeochron$x, xyGeochron$y, geochronMtEdgar$ANALYSIS_A) +
  coord_sf(expand =FALSE) 

# plot map with location of AMS data points
amsLocsOnKMap = pctKMap + geom_sf(data = amsData, aes(), size = 0.25, color = "white",  show.legend = FALSE) + coord_sf(expand = FALSE)  + annotate("rect", xmin = mtEdgarBboxUTM[[2]] - 21000, xmax = mtEdgarBboxUTM[[2]] - 1000, ymin = mtEdgarBboxUTM[[3]] + 1000, ymax = mtEdgarBboxUTM[[3]] + 4000, fill = "white") + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

ggsave("plotExports/amsLocationsMaps.pdf", plot = marrangeGrob(nrow = 1, ncol = 2, grobs = list(geochronOnGeoMap,amsLocsOnKMap), top = NULL), height = 4.2, unit = "in", useDingbats = FALSE)


# plot ams locations colored by unit
ggplot() + 
        geom_sf(data = amsDataMeans, aes(fill = SUPERSUITE, color = SUPERSUITE), shape = 21, size = 2) +
        scale_color_manual(values = superSuiteColors) +
        scale_fill_manual(values = superSuiteColors) +
        coord_sf(expand = FALSE, datum = crsUTM51S ) +
        theme_bw()


#labeled geologic map
unit_centers = st_geometry(st_centroid(mtEdgarMap[is.na(mtEdgarMap$SUPERSUITE) == FALSE,]))
xyLabel <- do.call(rbind, unit_centers) %>% data.frame() %>% setNames(c("x","y"))
xyLabel$unit = mtEdgarMap[is.na(mtEdgarMap$SUPERSUITE) == FALSE,]$UNITNAME
geoMapMtEdgar + geom_label(data = xyLabel, aes(x,y, label = unit), size = 2)




#-------------PLOT STRUCTURAL SYMBOLS ON MAP---------------------

# AMS FOLIATIONS AND LINEATIONS
dataSet = amsDataMeans

# Extract cartesian vectors of pole to foliation and direction of lineation
dataSet$pole = lapply(dataSet$rotation, function(rot) rot[3,])
dataSet$direction = lapply(dataSet$rotation, function(rot) rot[1,])

# Calculate dip and plunge values for the purpose of colorscale on map
dip = sapply(dataSet$pole, function(x) geoStrikeDipDegFromCartesian(x)[2])
plunge = sapply(dataSet$direction, function(x) geoTrendPlungeDegFromCartesian(x)[2])

# Plot geologic map with AMS location points colored by foliation dip (viridis color scale)
pSD_dip = geoMapMtEdgarSuperSuitesBW + 
  geom_sf(data = dataSet, aes(color = dip), size = 4) + scale_color_viridis_c(direction = -1) + 
  coord_sf(expand = FALSE)

# Plot strike and dip symbols on top of pSD_dip
pSD = plotMapSymbolsFromCartesian(plot = pSD_dip, dataSet, poleCartesian = dataSet$pole, symbolColor = "white") +  
  annotate("rect", xmin = mtEdgarBboxUTM[[2]] - 21000, xmax = mtEdgarBboxUTM[[2]] - 1000, ymin = mtEdgarBboxUTM[[3]] + 8000, ymax = mtEdgarBboxUTM[[3]] + 11000, fill = "black")

# Display plot
pSD 

# Plot geologic map with AMS location points colored by lineation plunge (viridis color scale)
pTP_plunge = geoMapMtEdgarSuperSuitesBW + geom_sf(data = dataSet, aes(color = plunge), size = 4) + scale_color_viridis_c(direction = -1) + coord_sf(expand = FALSE)

# Draw trend plunge symbols on the pTP_plunge map
pTP = plotMapSymbols(plot = pTP_plunge, dataSet, lineCartesian = dataSet$direction, symbolColor = "white") +
        annotate("rect", xmin = mtEdgarBboxUTM[[2]] - 21000, xmax = mtEdgarBboxUTM[[2]] - 1000, ymin = mtEdgarBboxUTM[[3]] + 8000, ymax = mtEdgarBboxUTM[[3]] + 11000, fill = "black") 

# Display plot
pTP

# Save maps to plotExports folder
ggsave("plotExports/mapStrikeDip.pdf", plot = pSD,  width = 5, height = 5, unit = "in")
ggsave("plotExports/mapTrendPlunge.pdf", plot = pTP,  width = 5, height = 5, unit = "in")



# FIELD FOLIATIONS
sdMap = plotMapSymbolsFromCartesian(plot = geoMapMtEdgarSuperSuitesBW, sf_dataFrame = fieldFabGranites2016, poleCartesian = fieldFabGranites2016$pole, symbolColor = "white")

sdMap 

sdMap2016_2018 = plotMapSymbols(plot = sdMap,sf_dataFrame = fieldFabGranites2018, poleCartesian = fieldFabGranites2018$pole, symbolColor = "white")

sdMap2016_2018 = sdMap2016_2018 + annotate("rect", xmin = mtEdgarBboxUTM[[2]] - 21000, xmax = mtEdgarBboxUTM[[2]] - 1000, ymin = mtEdgarBboxUTM[[3]] + 1000, ymax = mtEdgarBboxUTM[[3]] + 4000, fill = "black") + coord_sf(expand = FALSE, datum = crsDecDeg) 

ggsave("plotExports/fieldSdMap.pdf", plot = sdMap2016_2018, width = 4, height = 4, unit = "in")



#-----------PLOT MAPS OF AMS PARAMETERS----------------


# Shape parameter T
pT = formatParameterMapPlot(geoMapMtEdgarSuperSuitesBW + geom_sf(data = dataSet, aes(color = sapply(logA, ellLodeNu)), size = 4) + scale_color_gradient2(name = "T"))
pT
ggsave("plotExports/mapTparameter.pdf", plot = pT, height = 5,width = 5, unit = "in" )


# Mean susceptibility
pKm = formatParameterMapPlot(geoMapMtEdgarSuperSuitesBW + geom_sf(data = dataSet, aes(color = log10(Km)), size = 4) + scale_color_gradient2(name = "Km", midpoint = -3.3, low = muted("blue"), high = muted("red")))
pKm
ggsave("plotExports/mapKmparameter.pdf", plot = pKm, height = 5,width = 5, unit = "in" )


# Anisotropy degree Pj
pPj = formatParameterMapPlot(geoMapMtEdgarSuperSuitesBW + geom_sf(data = dataSet, aes(color = sapply(logA, ellJelinekP)), size = 4) + scale_color_viridis_c(name = "Pj", limits = c(1.0, 1.4), oob = squish))
pPj
ggsave("plotExports/mapPjparameter.pdf", plot = pPj, height = 5,width = 5, unit = "in" )



# -----------PLOT PRESENCE OF MAFICS/ACCESSORY PHASES-------------------------
microStructures$Biotite = grepl("Bi",microStructures$Mafics)
microStructures$Titanite = grepl("Ti",microStructures$Mafics)
microStructures$Hornblende = grepl("Hb",microStructures$Mafics)
microStructures$Oxides = grepl("Oxides",microStructures$Mafics)


Biotite = geoMapMtEdgarSuperSuitesBW + 
  geom_sf(data = microStructures, size = 0.5, color = "white") +
  geom_sf(data = microStructures[microStructures$Biotite == TRUE,], shape = 22, color = "white", fill = "black", size = 2) + 
  annotate("rect", xmin = mtEdgarBboxUTM[[2]] - 21000, xmax = mtEdgarBboxUTM[[2]] - 1000, ymin = mtEdgarBboxUTM[[3]] + 1000, ymax = mtEdgarBboxUTM[[3]] + 4000, fill = "black") + coord_sf(expand = FALSE) + theme_void() 

Hornblende = geoMapMtEdgarSuperSuitesBW + 
  geom_sf(data = microStructures, size = 0.5, color = "white") +
  geom_sf(data = microStructures[microStructures$Hornblende == TRUE,], shape = 22, color = "white", fill = "black", size = 2) +
  annotate("rect", xmin = mtEdgarBboxUTM[[2]] - 21000, xmax = mtEdgarBboxUTM[[2]] - 1000, ymin = mtEdgarBboxUTM[[3]] + 1000, ymax = mtEdgarBboxUTM[[3]] + 4000, fill = "black") + coord_sf(expand = FALSE) + theme_void()  

Titanite = geoMapMtEdgarSuperSuitesBW + 
  geom_sf(data = microStructures, size = 0.5, color = "white") +
  geom_sf(data = microStructures[microStructures$Titanite == TRUE,], shape = 22, color = "white", fill = "black", size = 2) + 
  annotate("rect", xmin = mtEdgarBboxUTM[[2]] - 21000, xmax = mtEdgarBboxUTM[[2]] - 1000, ymin = mtEdgarBboxUTM[[3]] + 1000, ymax = mtEdgarBboxUTM[[3]] + 4000, fill = "black") + coord_sf(expand = FALSE) + theme_void() 

Oxides = geoMapMtEdgarSuperSuitesBW + 
  geom_sf(data = microStructures, size = 0.5, color = "white") +
  geom_sf(data = microStructures[microStructures$Oxides == TRUE,], shape = 22, color = "white", fill = "black", size = 2) +
  annotate("rect", xmin = mtEdgarBboxUTM[[2]] - 21000, xmax = mtEdgarBboxUTM[[2]] - 1000, ymin = mtEdgarBboxUTM[[3]] + 1000, ymax = mtEdgarBboxUTM[[3]] + 4000, fill = "black") + coord_sf(expand = FALSE) + theme_void() 

ggsave("plotExports/biotiteMap.pdf", plot = Biotite, height = 3, width = 3, unit = "in", useDingbats = FALSE)
ggsave("plotExports/titaniteMap.pdf", plot = Titanite, height = 3, width = 3, unit = "in", useDingbats = FALSE)
ggsave("plotExports/hornblendeMap.pdf", plot = Hornblende, height = 3, width = 3, unit = "in", useDingbats = FALSE)
ggsave("plotExports/oxidesMap.pdf", plot = Oxides, height = 3, width = 3, unit = "in", useDingbats = FALSE)


#-------------PLOT MICROSTRUCTURAL CATEGORIZATIONSS-------------

microShapes = c("M" = 21,
                "SM2" = 22,
                "SM3" = 22,
                "HT" = 24)

microFill = c("M" = "black",
              "SM2" = "red",
              "SM3" = "pink",
              "HT" = "light blue")

microColor = c("M" = "white",
               "SM2" = "white",
               "SM3" = "black",
               "HT" = "black")

microStructures$Biotite = grepl("Bi",microStructures$Mafics)
microStructures$Titanite = grepl("Ti",microStructures$Mafics)
microStructures$Hornblende = grepl("Hb",microStructures$Mafics)
microStructures$Oxides = grepl("Oxides",microStructures$Mafics)


Biotite = geoMapMtEdgarSuperSuitesBW + 
  geom_sf(data = microStructures, size = 0.5, color = "white") +
  geom_sf(data = microStructures[microStructures$Biotite == TRUE,], shape = 22, color = "white", fill = "black", size = 2) + 
  annotate("rect", xmin = mtEdgarBboxUTM[[2]] - 21000, xmax = mtEdgarBboxUTM[[2]] - 1000, ymin = mtEdgarBboxUTM[[3]] + 1000, ymax = mtEdgarBboxUTM[[3]] + 4000, fill = "black") + coord_sf(expand = FALSE) + theme_void() 

Hornblende = geoMapMtEdgarSuperSuitesBW + 
  geom_sf(data = microStructures, size = 0.5, color = "white") +
  geom_sf(data = microStructures[microStructures$Hornblende == TRUE,], shape = 22, color = "white", fill = "black", size = 2) +
  annotate("rect", xmin = mtEdgarBboxUTM[[2]] - 21000, xmax = mtEdgarBboxUTM[[2]] - 1000, ymin = mtEdgarBboxUTM[[3]] + 1000, ymax = mtEdgarBboxUTM[[3]] + 4000, fill = "black") + coord_sf(expand = FALSE) + theme_void()  

Titanite = geoMapMtEdgarSuperSuitesBW + 
  geom_sf(data = microStructures, size = 0.5, color = "white") +
  geom_sf(data = microStructures[microStructures$Titanite == TRUE,], shape = 22, color = "white", fill = "black", size = 2) + 
  annotate("rect", xmin = mtEdgarBboxUTM[[2]] - 21000, xmax = mtEdgarBboxUTM[[2]] - 1000, ymin = mtEdgarBboxUTM[[3]] + 1000, ymax = mtEdgarBboxUTM[[3]] + 4000, fill = "black") + coord_sf(expand = FALSE) + theme_void() 

Oxides = geoMapMtEdgarSuperSuitesBW + 
  geom_sf(data = microStructures, size = 0.5, color = "white") +
  geom_sf(data = microStructures[microStructures$Oxides == TRUE,], shape = 22, color = "white", fill = "black", size = 2) +
  annotate("rect", xmin = mtEdgarBboxUTM[[2]] - 21000, xmax = mtEdgarBboxUTM[[2]] - 1000, ymin = mtEdgarBboxUTM[[3]] + 1000, ymax = mtEdgarBboxUTM[[3]] + 4000, fill = "black") + coord_sf(expand = FALSE) + theme_void() 

ggsave("plotExports/biotiteMap.pdf", plot = Biotite, height = 3, width = 3, unit = "in", useDingbats = FALSE)
ggsave("plotExports/titaniteMap.pdf", plot = Titanite, height = 3, width = 3, unit = "in", useDingbats = FALSE)
ggsave("plotExports/hornblendeMap.pdf", plot = Hornblende, height = 3, width = 3, unit = "in", useDingbats = FALSE)
ggsave("plotExports/oxidesMap.pdf", plot = Oxides, height = 3, width = 3, unit = "in", useDingbats = FALSE)


#--------------------------

microMap = geoMapMtEdgarSuperSuitesBW +
  geom_sf(data = amsDataMeans, size = 0.5, color = "white") +
  new_scale_fill() + new_scale_color() +
  geom_sf(data = microStructures, aes(shape = Classification, fill = Classification, color = Classification), size = 3, show.legend = "point") +
  scale_fill_manual(name = "Microstructure", values = microFill) +
  scale_shape_manual(name = "Microstructure",values = microShapes) +
  scale_color_manual(name = "Microstructure", values = microColor) + annotate("rect", xmin = mtEdgarBboxUTM[[1]] + 1000, xmax = mtEdgarBboxUTM[[1]] + 22000, ymin = mtEdgarBboxUTM[[3]] + 1000, ymax = mtEdgarBboxUTM[[3]] + 4000, fill = "white") + coord_sf(expand = FALSE, datum = crsDecDeg) +
  theme(panel.grid.major = element_blank(), 
        legend.position = c(1,0),
        legend.margin = margin(0.005,.005,.005,.005, "npc"),
        legend.justification = c(1,0),
        legend.key.height = unit(0.04, "npc"),
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white", color = "black", size = .5),
        legend.title = element_blank(),
        panel.border = element_rect(size = 1, color = "black", fill = NA))

microMap



ggsave("plotExports/microstructuresMap.pdf", plot = microMap, height = 4.5, width = 4.5, unit = "in", useDingbats = FALSE)




#-------PLOT MAP OF MAGNETIC MINERALOGY ANALYSIS LOCATIONS-------------------


magMinMap = geoMapMtEdgarSuperSuitesBW + 
  geom_sf(data = amsDataMeans, size = 2, color = "white", fill = "White") +
  geom_sf(data = magMinLocations[magMinLocations$VSM == 1,], shape = 24, color = "black", size = 1) +
  geom_sf(data = magMinLocations[magMinLocations$ThermSusc == 1,], shape = 21, color = "black", size = 2) +
  annotate("rect", xmin = mtEdgarBboxUTM[[2]] - 21000, xmax = mtEdgarBboxUTM[[2]] - 1000, ymin = mtEdgarBboxUTM[[3]] + 1000, ymax = mtEdgarBboxUTM[[3]] + 4000, fill = "black") + coord_sf(expand = FALSE, datum = crsDecDeg) 

ggsave("plotExports/magMinMap.pdf", plot = magMinMap, height = 4.5, width = 4.5, unit = "in", useDingbats = FALSE)

#--------PLOT MAP OF TRIAXIAL, PROLATE-GIRDLED, AND OBLATE-GIRDLED STATIONSS ------------------



triaxialityShape = c("Triaxial" = 23,
                     "Prolate" = 25,
                     "Oblate" = 21,
                     "REJECT" = 4)

triaxialityFill = c("Triaxial" = "black",
                    "Prolate" = "red",
                    "Oblate" = "blue",
                    "REJECT" = "grey")

geoMapMtEdgarSuperSuitesBW +
  new_scale_fill()+
  geom_sf(data = ellipsoidTriaxiality, aes(fill = Triaxiality), color = "white", shape = 21, show.legend = "point", size = 4) +
  scale_fill_manual(values = triaxialityFill)

        