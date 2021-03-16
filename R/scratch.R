

# 
# plotEqualAreaNet(points = list(lapply(amsDataMeans$rotation, function(s) s[1,])),
#                  colorBy = sapply(amsDataMeans$logA, ellLodeNu), showLegend = TRUE) + scale_fill_viridis_c()
# 
# plotEqualAreaNet(points = list(lapply(amsDataMeans$rotation, function(s) s[3,])),
#                  colorBy = sapply(amsDataMeans$logA, ellLodeNu), showLegend = TRUE, shape = 21) + scale_fill_viridis_c()
# 
# 
# plotEqualAreaNet(points = list(lapply(amsDataMeans$rotation, function(s) s[1,])),
#                  colorBy = amsDataMeans$Triaxiality, showLegend = TRUE) + scale_fill_manual(values = triaxialityFill)
# 
# plotEqualAreaNet(points = list(lapply(amsDataMeans$rotation, function(s) s[3,])),
#                  colorBy = amsDataMeans$Triaxiality, showLegend = TRUE) + scale_fill_manual(values = triaxialityFill)

















compPlots = list()
count = 1
maxLim = sqrt(max(unlist(amsDataMeans$vector)^2))

for(i in 1:5){
    for(n in 1:5){
        ellComponents = data.frame(cbind(sapply(amsDataMeans$vector, function(s) s[n]),
                                         sapply(amsDataMeans$vector, function(s) s[i])))
        compPlots[[count]] = ggplot() + 
            geom_point(data = ellComponents, aes(x= X1, y = X2)) +
            coord_fixed(ratio = 1, xlim = c(-maxLim, maxLim), ylim = c(-maxLim, maxLim)) +
            theme_bw()+
            xlab(paste("comp.",n)) +
            ylab(paste("comp.",i)) +
            theme(panel.grid = element_blank())
        # axis.title = element_blank()) 
        count = count + 1
    }
}

plot_grid(plotlist = compPlots)




grid.arrange(grobs = compPlots, nrow = 5)



MeanEllChar = read.csv("data/MeanEllipsoidCategorization.csv")





















vectors = t(sapply(amsDataMeans$vector, function(s) s))
clusters = kmeans(vectors,5, iter.max = 100000)


withinClusterSS = c()
for (i in 2:20) {
    cluster = kmeans(vectors, i, iter.max =100000)
    withinClusterSS[i-1] = mean(cluster$withinss)
}

withinClustersSSFrame = data.frame("number" = c(2:20), "withinClusterSS" = withinClusterSS)

ggplot(data = withinClustersSSFrame, aes(x = number, y = withinClusterSS)) +
    geom_point() 

clusters = kmeans(vectors,5, iter.max = 20)
amsDataMeans$cluster = clusters$cluster
amsDataMeans$clusterOctNorm = clustersNorm$cluster


geoMapMtEdgarSuperSuitesBW +
    new_scale_fill() +
    geom_sf(data =amsDataMeans, aes(fill = factor(cluster)), size = 2,shape = 21, show.legend = "point") +
    theme_bw()



plotEqualAreaNet(points = list(lapply(amsDataMeans$rotation, function(s) s[1,])),
                 colorBy = factor(clusters$cluster), showLegend = TRUE)


plotEqualAreaNet(points = list(lapply(amsDataMeans$rotation, function(s) s[3,])),
                 colorBy = factor(clusters$cluster))


plotPjT() + geom_point(data = amsDataMeans, aes(x = sapply(logA, ellJelinekP), y = sapply(logA, ellLodeNu), color = factor(cluster)),shape = 19,size = 3, alpha = 1)



geoMapMtEdgarSuperSuitesBW +
    new_scale_fill() 
    #geom_sf(data =amsDataMeans, aes(fill = factor(cluster)), size = 2,shape = 21, show.legend = "point") +
    ggplotly(ggplot() + geom_text(data= amsDataMeans, aes(x = easting, y = northing, label = station.x )) +
    theme_bw())








K1 = function(x) {x[1,]}
K3 = function(x) {x[3,]}



marginAMS = amsData[amsData$station %in% list("AME18_017", "AME18_104", "AME16_141", "AME18_015", "AME18_242") == TRUE,]

plotEqualAreaNet(points = list(lapply(marginAMS$rotation, K1)),
                 curves = list(lapply(marginAMS$rotation, K3)))



amsData[amsData$station == "AME18_018",]













dataSet = amsDataMeans





# strike dip trend plunge for AMS data
sd <- lapply(dataSet$rotation, function(s) geoStrikeDipDegFromCartesian(s[3,]))
tp <- lapply(dataSet$rotation, function(s) geoTrendPlungeDegFromCartesian(lower(s[1,])))

plot = ggplot() +
    xlim(c(mtEdgarBboxUTM[1], mtEdgarBboxUTM[2])) +
    ylim(c(mtEdgarBboxUTM[3], mtEdgarBboxUTM[4]))

xyCoords <- do.call(rbind, st_geometry(dataSet)) %>% 
    as_tibble() %>% setNames(c("x","y","elev"))




# define strike dip symbol
# sdLine = data.frame(x = c(0,0),
#                     y = c(-1,1))
# sdTriangle = data.frame(x = c(0,0.25,0,0),
#                         y = c(0.25, 0, -0.25, 0.25))
sdTriangle = data.frame(x = c(0,0,0,.25,0,0),
                        y = c(-1,1,0.25,0,-.25,-1))
ggplot() +
    #geom_line(data = sdLine,aes(x = x, y =y)) +
    geom_polygon(data = sdTriangle, aes(x=x, y=y), color = "black", fill = "black") +
    coord_fixed(ratio = 1)

tpPolygon = data.frame(x =c(0,0,.2,-.2,0,0), 
                       y =c(-1,1.4,1,1,1.4,0))
ggplot() +
    #geom_line(data = sdLine,aes(x = x, y =y)) +
    geom_polygon(data = tpPolygon, aes(x=x, y=y), color = "black", fill = "black") +
    coord_fixed(ratio = 1)


#Draw the plot
pSD = pctKMap + 
    geom_sf(data = dataSet, aes(color = dip), size = 4) + scale_color_viridis_c(direction = -1) + 
    coord_sf(expand = FALSE)



symbolColor = "black"

#rotate the symbols and plot them, one by one
for(i in 1:length(sd)){
    symbol = rotateSymbol(sdTriangle, sd[[i]][1], halfSize = 1000)
    # symbol[[1]]$X1 = symbol[[1]]$X1 + xyCoords$x[i]
    # symbol[[1]]$X2 = symbol[[1]]$X2 + xyCoords$y[i]
    symbol$X1 = symbol$X1 + xyCoords$x[i]
    symbol$X2 = symbol$X2 + xyCoords$y[i]
    pSD = pSD + 
        #geom_line(data = symbol[[1]],aes(x = X1, y =X2), color = symbolColor) +
        geom_polygon(data = symbol, aes(x=X1, y=X2), color = symbolColor, fill = symbolColor, size = .5)
}
pSD = pSD + theme_void() + 
    theme(panel.grid.major = element_blank(), 
          legend.position = c(1,0),
          legend.margin = margin(0.005,.005,.005,.005, "npc"),
          legend.justification = c(1,0),
          legend.key.height = unit(0.04, "npc"), 
          legend.direction = "horizontal",
          legend.background = element_rect(fill = "white", color = "black", size = .5),
          panel.border = element_rect(size = 1, color = "black", fill = NA))

pSD

#Draw the plot
pTP = pctKMap + geom_sf(data = sfDataF, aes(color = plunge), size = 4) + scale_color_viridis_c(direction = -1) + coord_sf(expand = FALSE)
symbolColor = "white"

#rotate the symbols and plot them, one by one
for(i in 1:length(tp)){
    
    arrowAdj = (tp[[i]][2]/90)
    tpPolygon = data.frame(x =c(0,0,.2,-.2,0,0), 
                           y =c(-1 + arrowAdj,1.4-arrowAdj,1-arrowAdj,1 -arrowAdj,1.4 - arrowAdj,-1 + arrowAdj))
    symbol = rotateSymbol(tpPolygon, tp[[i]][1], halfSize = 1000)
    # symbol[[1]]$X1 = symbol[[1]]$X1 + xyCoords$x[i]
    # symbol[[1]]$X2 = symbol[[1]]$X2 + xyCoords$y[i]
    symbol$X1 = symbol$X1 + xyCoords$x[i]
    symbol$X2 = symbol$X2 + xyCoords$y[i]
    pTP = pTP + 
        #geom_line(data = symbol[[1]],aes(x = X1, y =X2), color = symbolColor) +
        geom_polygon(data = symbol, aes(x=X1, y=X2), color = symbolColor, fill = symbolColor)
}

pTP = pTP + theme_void() + 
    theme(panel.grid.major = element_blank(), 
          legend.position = c(1,0),
          legend.margin = margin(0.005,.005,.005,.005, "npc"),
          legend.justification = c(1,0),
          legend.key.height = unit(0.04, "npc"), 
          legend.direction = "horizontal",
          legend.background = element_rect(fill = "white", color = "black", size = .5),
          panel.border = element_rect(size = 1, color = "black", fill = NA))

ggsave("plotExports/mapStrikeDip.pdf", plot = pSD,  width = 5, height = 5, unit = "in")
ggsave("plotExports/mapTrendPlunge.pdf", plot = pTP,  width = 5, height = 5, unit = "in")



#Draw the plot
pT = geoMapMtEdgarSuperSuitesBW + geom_sf(data = dataSet, aes(color = sapply(logA, ellLodeNu)), size = 4) + scale_color_gradient2(name = "T") + coord_sf(expand = FALSE) + theme_void() + 
    theme(panel.grid.major = element_blank(), 
          legend.position = c(1,0),
          legend.margin = margin(0.005,.005,.005,.005, "npc"),
          legend.justification = c(1,0),
          legend.key.height = unit(0.04, "npc"), 
          legend.direction = "horizontal",
          legend.background = element_rect(fill = "white", color = "black", size = .5),
          panel.border = element_rect(size = 1, color = "black", fill = NA)) 

ggsave("plotExports/mapTparameter.pdf", plot = pT, height = 5,width = 5, unit = "in" )








#comparison betwen AMS and field fabric

stationFolComp = data.frame(Name = factor(),Distance = double(),PValue = double())

for(i in 1:length(fieldFabCores$Name)){
    station = fieldFabCores[i,]
    stationAMS = amsData[amsData$station == station$Name,]
    
    if(is.na(stationAMS$Name[1]) == FALSE){
        
        boots = lineBootstrapInference(lapply(stationAMS$rotation, function(s) s[3,]), 10000)
        distance = lineDistance(geoCartesianFromStrikeDipDeg(c(station$strike[1], station$dip[1])), boots$center)/degree
        
        p = boots$pvalue(geoCartesianFromStrikeDipDeg(c(station$strike[1], station$dip[1])))
        Name = station$Name
        newRow = data.frame(Name = station$Name, Distance = distance,PValue = p)
        stationFolComp = rbind(stationFolComp, newRow)
        
        # stationAMSPlot = plotEqualAreaNet(points = list(boots$us, lapply(stationAMS$rotation, function(s) s[3,])), pointSize = list(1,3))
        # 
        # plot = plotEqualAreaNet(plot = stationAMSPlot, points = list(list(replicate(2,geoCartesianFromStrikeDipDeg(c(station$strike, station$dip))))), color = "orange")
    }
    # print(plot)
}

fieldFab2018AMSstations = fieldFabGranites2018[is.na(fieldFabGranites2018$sampleName) == FALSE,]

fieldFab2018AMSstations$Name = fieldFab2018AMSstations$sampleName

for(i in 1:length(fieldFab2018AMSstations$Name)){
    station = fieldFab2018AMSstations[i,]
    stationAMS = amsData[amsData$station == substr(as.character(station$Name), 1, 9),]
    
    if(is.na(stationAMS$Name[1]) == FALSE){
        
        boots = lineBootstrapInference(lapply(stationAMS$rotation, function(s) s[3,]), 10000)
        distance = lineDistance(geoCartesianFromStrikeDipDeg(c(station$strike[1], station$dip[1])), boots$center)/degree
        
        p = boots$pvalue(geoCartesianFromStrikeDipDeg(c(station$strike[1], station$dip[1])))
        Name = station$Name
        newRow = data.frame(Name = substr(as.character(station$Name), 1, 9), Distance = distance,PValue = p)
        stationFolComp = rbind(stationFolComp, newRow)
        
         stationAMSPlot = plotEqualAreaNet(points = list(boots$us, lapply(stationAMS$rotation, function(s) s[3,])), pointSize = list(1,3))
         
         plot = plotEqualAreaNet(plot = stationAMSPlot, points = list(list(replicate(2,geoCartesianFromStrikeDipDeg(c(station$strike, station$dip))))), color = "orange")
    }
     print(plot)
}

ggplot() +
    geom_histogram(data = stationFolComp, aes(x = Distance))

ggplot() +
    geom_histogram(data = stationFolComp, aes(x = PValue))





stationFolComp
stationFolComp$Station = stationFolComp$Name

stationFolComp = locationPairing(stationFolComp, locData = locData, 9 )




