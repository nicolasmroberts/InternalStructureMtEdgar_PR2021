
# JOORINA ARCH WORK

joorinaArch = amsDataMeans[amsDataMeans$UNITNAME == "Joorina Granodiorite" | amsDataMeans$UNITNAME == "Kennell Granodiorite" | amsDataMeans$UNITNAME == "Bishop Creek Monzogranite",]

sd = lapply(1:length(joorinaArch$station.x), function(s) geoStrikeDipDegFromCartesian(joorinaArch$rotation[[s]][3,]))
dip = sapply(sd, function(s) s[2])

tp = lapply(1:length(joorinaArch$station.x), function(s) geoTrendPlungeDegFromCartesian(lower(joorinaArch$rotation[[s]][1,])))
plunge = sapply(tp, function(s) s[2])


JoorinaSD = plotMapSymbols(plot = geoMapMtEdgarSuperSuitesBW, sf_dataFrame = joorinaArch, poleCartesian = lapply(joorinaArch$rotation, function(s) s[3,]), colorBy = dip, pointSize = 10) + coord_sf(xlim = c(175000, 200000), ylim = c(7635000, 7660000))


JoorinaTP = plotMapSymbols(plot = geoMapMtEdgarSuperSuitesBW, sf_dataFrame = joorinaArch, lineCartesian = lapply(joorinaArch$rotation, function(s) s[1,]), colorBy = plunge, pointSize = 10) + 
     coord_sf(xlim = c(175000, 200000), ylim = c(7635000, 7660000), datum = crsUTM51S)  
     #coord_sf(xlim = c(190000, 200000), ylim = c(7650000, 7660000), datum = crsUTM51S)  
      #geom_sf_label(data = joorinaArch, aes(label = station.x))


JoorinaSD
JoorinaTP

ggsave("plotExports/JoorinaSD.pdf", plot = JoorinaSD, height = 4.2, width = 4.2, unit = "in")
ggsave("plotExports/JoorinaTP.pdf", plot = JoorinaTP, height = 4.2, width = 4.2, unit = "in")


JoorinaTip = amsData[amsData$station == "AME16_144" | amsData$station == "AME16_142"|amsData$station == "AME18_018",]

JoorinaMiddle = amsData[amsData$station == "AME16_090" | amsData$station == "AME16_091" | amsData$station == "AME18_165" | amsData$station == "AME18_167" | amsData$station == "AME18_168" | amsData$station == "AME18_169",]

JoorinaTail = amsData[amsData$station == "AME16_079" | amsData$station == "AME18_130"|amsData$station == "AME18_131" | amsData$station == "AME16_080" | amsData$station == "AME18_128",] 



plotTip = plotEqualAreaNet(points = list(lapply(JoorinaTip$rotation, function(s) s[1,])))
plotMiddle = plotEqualAreaNet(points = list(lapply(JoorinaMiddle$rotation, function(s) s[1,])))
plotTail = plotEqualAreaNet(points = list(lapply(JoorinaTail$rotation, function(s) s[1,])))

ggsave("plotExports/plotTip.pdf", plot = plotTip, width = 2.5, height = 2.5, unit = "in")
ggsave("plotExports/plotMiddle.pdf", plot = plotMiddle, width = 2.5, height = 2.5, unit = "in")
ggsave("plotExports/plotTail.pdf", plot = plotTail, width = 2.5, height = 2.5, unit = "in")




# TAMBINA LOBE WORK

sd = lapply(1:length(amsDataMeans$station.x), function(s) geoStrikeDipDegFromCartesian(amsDataMeans$rotation[[s]][3,]))
dip = sapply(sd, function(s) s[2])

tp = lapply(1:length(amsDataMeans$station.x), function(s) geoTrendPlungeDegFromCartesian(lower(amsDataMeans$rotation[[s]][1,])))
plunge = sapply(tp, function(s) s[2])

WLobeSD = plotMapSymbols(plot = geoMapMtEdgarSuperSuitesBW, sf_dataFrame = amsDataMeans, poleCartesian = lapply(amsDataMeans$rotation, function(s) s[3,]), colorBy = dip, pointSize = 12) + 
     coord_sf(xlim = c(172000, 192000), ylim = c(7645000,7670000))

WLobeTP = plotMapSymbols(plot = geoMapMtEdgarSuperSuitesBW, sf_dataFrame = amsDataMeans, lineCartesian = lapply(amsDataMeans$rotation, function(s) s[1,]), colorBy = plunge, pointSize = 12) + 
     coord_sf(xlim = c(172000, 192000), ylim = c(7645000,7670000)) 
     # geom_sf_label(data = amsDataMeans, aes(label = station.x))


WLobeSD
WLobeTP

ggsave("plotExports/WLobSD.pdf", plot = WLobeSD, width = 3.8, height = 4.6, unit = "in")
ggsave("plotExports/WLobTP.pdf", plot = WLobeTP, width = 3.8, height = 4.6, unit = "in")


TamWLobe1 = amsData[amsData$station == "AME16_145" | amsData$station == "AME16_146"  | amsData$station == "AME16_089" | amsData$station == "AME16_078" | amsData$station == "AME16_077" | amsData$station == "AME16_076" | amsData$station == "AME18_265" | amsData$station == "AME16_145" | amsData$station == "AME16_075",]

TamWLobe2 = amsData[amsData$station == "AME16_161" | amsData$station == "AME16_162"  | amsData$station == "AME16_056" | amsData$station == "AME16_055" | amsData$station == "AME16_035" | amsData$station == "AME16_034" ,]

TamWLobeHinge = amsData[amsData$station == "AME16_029" | amsData$station == "AME16_030"  | amsData$station == "AME16_031" | amsData$station == "AME16_032" | amsData$station == "AME16_033" ,]

plotTWL1 = plotEqualAreaNet(points = list(lapply(TamWLobe1$rotation, function(s) s[3,])), shape = 21)
plotHinge = plotEqualAreaNet(points = list(lapply(TamWLobeHinge$rotation, function(s) s[3,])), shape = 21)
plotTWL2 = plotEqualAreaNet(points = list(lapply(TamWLobe2$rotation, function(s) s[3,])), shape = 21)


plotTWL1
plotHinge
plotTWL2


ggsave("plotExports/plotTWL1.pdf", plot = plotTWL1, width = 2.5, height = 2.5, unit = "in")
ggsave("plotExports/plotHinge.pdf", plot = plotHinge, width = 2.5, height = 2.5, unit = "in")
ggsave("plotExports/plotTWL2.pdf", plot = plotTWL2, width = 2.5, height = 2.5, unit = "in")




SWLobeSD = plotMapSymbols(plot = geoMapMtEdgarSuperSuitesBW, sf_dataFrame = amsDataMeans, poleCartesian = lapply(amsDataMeans$rotation, function(s) s[3,]), colorBy = dip, pointSize = 12) + 
     coord_sf(xlim = c(172000, 198000), ylim = c(7630000,7645000))
SWLobeTP = plotMapSymbols(plot = geoMapMtEdgarSuperSuitesBW, sf_dataFrame = amsDataMeans, lineCartesian = lapply(amsDataMeans$rotation, function(s) s[1,]), colorBy = plunge, pointSize = 12) + 
      coord_sf(xlim = c(172000, 198000), ylim = c(7630000,7645000))  
# geom_sf_label(data = amsDataMeans, aes(label = station.x))

ggsave("plotExports/SWLobSD.pdf", plot = SWLobeSD, width = 5, height = 3, unit = "in")
ggsave("plotExports/SWLobTP.pdf", plot = SWLobeTP, width = 5, height = 3, unit = "in")




TamSWLobeS = amsData[amsData$station == "AME16_137" | amsData$station == "AME18_242" | amsData$station == "AME16_138" | amsData$station == "AME18_111",]
TamSWLobeN = amsData[amsData$station == "AME16_139" | amsData$station == "AME16_140" | amsData$station == "AME18_114",]

plotTSWLobeS = plotEqualAreaNet(points = list(lapply(TamSWLobeS$rotation, function(s) s[3,])), shape = 21)
plotTSWLobeN = plotEqualAreaNet(points = list(lapply(TamSWLobeN$rotation, function(s) s[3,])), shape = 21)

plotTSWLobeN
plotTSWLobeS

ggsave("plotExports/plotTSWLobeN.pdf", plot = plotTSWLobeN, width = 2.5, height = 2.5, unit = "in")
ggsave("plotExports/plotTSWLobeS.pdf", plot = plotTSWLobeS, width = 2.5, height = 2.5, unit = "in")

# CLELAND WORK
clelandTP = plotMapSymbols(plot = geoMapMtEdgarSuperSuitesBW, sf_dataFrame = amsDataMeans, lineCartesian = lapply(amsDataMeans$rotation, function(s) s[1,]), colorBy = dip, pointSize = 7, symbolColor = "white") + 
     coord_sf(xlim = c(185000, 220000), ylim = c(7645000,7675000))
clelandTP

ggsave("plotExports/celandTP.pdf", plot = clelandTP, width = 4.5, height = 4, unit = "in")



