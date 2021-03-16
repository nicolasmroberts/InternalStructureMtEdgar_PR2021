stationList <- amsData$station
# Remove duplicate station names
stationList <- unique(stationList, incomparables = FALSE)

# Automatically plot and save equal area nets and Pj v. T plots for all stations
count = 0
for(i in 1:length(stationList)){
     stationName = stationList[i]
     plotEqualAreaPjT(amsData, stationName, showBoots = TRUE, showCI = FALSE, numBoots = 1000, savePDF = TRUE, saveFolderPath = "plotExports/")
     count = count + 1
     print(count)
}


# Automatically plot and save Pj v. T plots for each stations
count = 0
for(i in 1:length(stationList)){
     stationName = stationList[[i]]
     bootTest = computeAMSBootStrap(amsData[amsData$station == stationName,], numBoots = 1000)
     if(class(bootTest) != "character"){
          KmvPj = plotKmvPj(amsData = amsData[amsData$station == stationName,], colorBySupersuite = FALSE, bootsFrame = bootTest,  colorKey = superSuiteColors)
     } else {KmvPj = blankplotKmvPj()}
     KmvPj = plotKmvPj(plot=KmvPj, amsData = amsDataMeans[amsDataMeans$station == stationName,], colorBySupersuite = FALSE, bootsFrame = NULL,  colorKey = superSuiteColors, pointSize = 3) + ylim(c(1,3)) + scale_x_log10(limits = c(1e-6,1e-2))
     ggsave(paste0("plotExports/",stationName,"_KmvPj.pdf"), plot = KmvPj, width = 4, height = 2.5, unit = "in")
     count = count + 1
     print(count)
}


# set a station for individual plots
stationName = "AME18_131"

# plot equal area net and Pj v. T for an individual station
plotEqualAreaPjT(amsData, stationName, showBoots = TRUE, showCI = FALSE, numBoots = 1000, savePDF = FALSE)

# plot Km v Pj for an individual station
bootTest = computeAMSBootStrap(amsData[amsData$station == stationName,], numBoots = 1000)
KmvPj = plotKmvPj(amsData = amsData[amsData$station == stationName,], colorBySupersuite = FALSE, bootsFrame = bootTest,  colorKey = superSuiteColors) + ylim(c(1,2)) + scale_x_log10(limits = c(1e-6,1e-2))
KmvPj = plotKmvPj(plot=KmvPj, amsData = amsDataMeans[amsDataMeans$station == stationName,], colorBySupersuite = FALSE, bootsFrame = NULL,  colorKey = superSuiteColors, pointSize = 3)
KmvPj





# Automatically plot and save equal area nets and Pj v. T plots for all stations
count = 0
plots = list()
for(i in 1:length(stationList)){
        stationName = stationList[i]
        
        t = ggplot() + annotate("text", x = 0, y = 0, angle = 90, label = stationName) + theme_void()
        p = plotEqualAreaPjT(amsData, stationName, showBoots = TRUE, showCI = FALSE, numBoots = 1000, savePDF = FALSE)
        plots = c(plots, list(t), p)
        count = count + 1
        print(count)
}



layMatrix = rbind(c(1,2,2,3,3,3),
                  c(4,5,5,6,6,6),
                  c(7,8,8,9,9,9),
                  c(10,11,11,12,12,12),
                  c(13,14,14,15,15,15))
                  

ggsave("plotExports/EANandPjTPlots.pdf", plot = marrangeGrob(grobs = plots, layout_matrix = layMatrix, top = NULL, heights = unit(c(1.75,1.75,1.75,1.75,1.75), c("in","in","in","in","in")), widths = unit(c(.5,1,1,1,1,1), c("in","in","in","in","in","in"))), width = 8.5, height = 11, unit = "in")





