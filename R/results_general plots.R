# dataAMS = amsData
dataAMS = amsDataMeans
# ------------------------------------------
# A series of common 2D plots for all the data

plotHsu(dataAMS, plotMean = FALSE)
plotShapeVersusAnisotropy(dataAMS, gradient = "none")

plotShapeVersusKm(dataAMS,gradient = "none")
plotOctVersusKm(dataAMS, gradient = "none")

plotHistogramLodes(dataAMS, binWidth = 0.1)
plotHistogramOct(dataAMS, binWidth = 0.02)


# --------------PARAMETER PLOTS (ALL STATIONS)----------

# Susceptibility vs. Anisotropy Pj
allCoresKmvPj  = plotKmvPj(amsData=amsData, pointSize = 0.5, pointFill = "grey", pointEdge = "grey") + geom_vline(xintercept = 10^-3.31)
meanKmvPj = plotKmvPj(plot=allCoresKmvPj, amsData = amsDataMeans, pointSize = 0.5, pointFill = "black" ) + theme(legend.position = "none")

# Anisotropy Pj vs. shape parameter T
allCoresPjT = plotPjT(amsData = amsData, pointSize = 0.5, pointFill = "grey", pointEdge = "grey")
meanPjT =plotPjT(plot = allCoresPjT, amsData= amsDataMeans, pointSize = 0.5, pointFill = "black") + theme(legend.position = "none") + xlim(c(1,2))

# Susceptibility vs. shape parameter T
allCoresKmT = plotKmT(amsData= amsData, pointSize = 0.5, pointFill = "grey", pointEdge = "grey")
meanKmT = plotKmT(plot = allCoresKmT, amsData = amsDataMeans, pointSize = 0.5, pointFill = "black") + theme(legend.position = "none")

ggsave("plotExports/KmvPjAllCores.pdf", plot = allCoresKmvPj, height = 2.5, width = 2.5, unit = "in", useDingbats = FALSE)
ggsave("plotExports/PjTAllCores.pdf", plot = allCoresPjT, height = 2.5, width = 2.5, unit = "in", useDingbats = FALSE)
ggsave("plotExports/KmTAllCores.pdf", plot = allCoresKmT, height = 2.5, width = 2.5, unit = "in", useDingbats = FALSE)
ggsave("plotExports/meanPjT.pdf", plot = meanPjT,  height = 2.5, width = 2.5, unit = "in", useDingbats = FALSE)
ggsave("plotExports/meanKmvPj.pdf", plot = meanKmvPj,  height = 2.5, width = 2.5, unit = "in", useDingbats = FALSE)
ggsave("plotExports/meanKmT.pdf", plot = meanKmT,  height = 2.5, width = 2.5, unit = "in", useDingbats = FALSE)



# Kernel density plots of different parameters. 

amsData[is.na(amsData$SUPERSUITE) == TRUE,]$SUPERSUITE = "Callina Supersuite" 
amsDataMeans[is.na(amsDataMeans$SUPERSUITE) == TRUE,]$SUPERSUITE = "Callina Supersuite" 

shapeParameterBySupersuite = ggplot() +
    geom_density(kernel = "gaussian", data = amsDataMeans[amsDataMeans$SUPERSUITE == "Tambina Supersuite"|amsDataMeans$SUPERSUITE == "Emu Pool Supersuite"| amsDataMeans$SUPERSUITE == "Cleland Supersuite",], aes(x= sapply(logA, ellLodeNu), fill = SUPERSUITE), alpha = .75, show.legend = FALSE)+
    scale_fill_manual(values = superSuiteColorsBW) +
    coord_cartesian(expand=FALSE) +
    xlab("T") +
    ylab("Density") +
    ylim(c(0,1)) + 
    geom_vline(xintercept = 0) +
    theme_classic()

PjParameterBySupersuite = ggplot() +
    geom_density(kernel = "gaussian", data = amsDataMeans[amsDataMeans$SUPERSUITE == "Tambina Supersuite" | amsDataMeans$SUPERSUITE == "Emu Pool Supersuite"| amsDataMeans$SUPERSUITE == "Cleland Supersuite",], aes(x= sapply(logA, ellJelinekP), fill = SUPERSUITE), alpha = .75, show.legend = FALSE)+
    scale_fill_manual(values = superSuiteColorsBW) +
    coord_cartesian(expand=FALSE) +
    xlab("Pj") +
    ylab("Density") +
    xlim(c(1,2)) +
    theme_classic()

KmParameterBySupersuite = ggplot() +
    geom_density(kernel = "gaussian", data = amsDataMeans[ amsDataMeans$SUPERSUITE == "Tambina Supersuite"|amsDataMeans$SUPERSUITE == "Emu Pool Supersuite"| amsDataMeans$SUPERSUITE == "Cleland Supersuite",], aes(x=Km, fill = SUPERSUITE), alpha = .75, show.legend = FALSE)+
    scale_fill_manual(values = superSuiteColorsBW) +
    coord_cartesian(expand=FALSE) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)))  +  
    xlab("Bulk Susceptibility (SI)") +        geom_vline(xintercept = 0) +
    ylab("Density") +
    ylim(c(0,1)) + 
    geom_vline(xintercept = 0) +
    theme_classic()


# Kernel density plots of different parameters. 
shapeParameterBySupersuiteALL = ggplot() +
    geom_density(kernel = "gaussian", data = amsData[amsData$SUPERSUITE == "Tambina Supersuite"|amsData$SUPERSUITE == "Emu Pool Supersuite"| amsData$SUPERSUITE == "Cleland Supersuite",], aes(x= sapply(logA, ellLodeNu), fill = SUPERSUITE), alpha = .75, show.legend = FALSE)+
    scale_fill_manual(values = superSuiteColorsBW) +
    coord_cartesian(expand=FALSE) +
    xlab("T") +
    ylab("Density") +
    ylim(c(0,1)) + 
    geom_vline(xintercept = 0) +
    theme_classic()

PjParameterBySupersuiteALL = ggplot() +
    geom_density(kernel = "gaussian", data = amsData[amsData$SUPERSUITE == "Tambina Supersuite" | amsData$SUPERSUITE == "Emu Pool Supersuite"| amsData$SUPERSUITE == "Cleland Supersuite",], aes(x= sapply(logA, ellJelinekP), fill = SUPERSUITE), alpha = .75,show.legend = FALSE)+
    scale_fill_manual(values = superSuiteColorsBW) +
    coord_cartesian(expand=FALSE) +
    xlab("Pj") +
    ylab("Density") +
    xlim(c(1,2)) +
    theme_classic()

KmParameterBySupersuiteALL = ggplot() +
    geom_density(kernel = "gaussian", data = amsData[ amsData$SUPERSUITE == "Tambina Supersuite"|amsData$SUPERSUITE == "Emu Pool Supersuite"| amsData$SUPERSUITE == "Cleland Supersuite",], aes(x=Km, fill = SUPERSUITE), alpha = .75,show.legend = FALSE)+
    scale_fill_manual(values = superSuiteColorsBW) +
    coord_cartesian(expand=FALSE) +
    ylab("Density") +
    ylim(c(0,1)) + 
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)))  +  
    xlab("Bulk Susceptibility (SI)") +        geom_vline(xintercept = 0) +
    theme_classic()


ggsave("plotExports/ShapeParameterDensityGraph.pdf", plot = shapeParameterBySupersuite, height = 2.5, width = 2.5, unit = "in")
ggsave("plotExports/KmParameterDensityGraph.pdf", plot = KmParameterBySupersuite, height = 2.5, width = 2.5, unit = "in")
ggsave("plotExports/PjParameterDensityGraph.pdf", plot = PjParameterBySupersuite, height = 2.5, width = 2.5, unit = "in")
ggsave("plotExports/ShapeParameterDensityGraphALL.pdf", plot = shapeParameterBySupersuiteALL, height = 3, width = 3.5, unit = "in")
ggsave("plotExports/KmParameterDensityGraphALL.pdf", plot = KmParameterBySupersuiteALL, height = 3, width = 3.5, unit = "in")
ggsave("plotExports/PjParameterDensityGraphALL.pdf", plot = PjParameterBySupersuiteALL, height = 3, width = 3.5, unit = "in")



# Susceptibility histogram for all cores (grey) and station mean (black)
meanKm = data.frame(Km = ((unlist(unlist(lapply(amsDataMeans$Km, function(s) c(s,s,s)))))))

susceptHist = ggplot() +
    geom_histogram(data = amsData, aes(x=Km), show.legend = FALSE, color = "white", fill="grey80")+
    scale_y_continuous(sec.axis = sec_axis(~ . /3)) +
    geom_histogram(data =meanKm, aes(x=Km), alpha = 1, color = "white", fill = "black") +
    coord_cartesian(expand=FALSE) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)))  +  
    xlab("Bulk Susceptibility (SI)") +  
    ylab(label = element_blank()) +
    theme_classic()

ggsave("plotExports/susceptHist.pdf", plot = susceptHist,  height = 3, width = 4, unit = "in", useDingbats = FALSE)





# ----------FIELD FOLIATION EQUAL AREA NETS--------------------------------

# Plot all field foliations together
# Combine data from two field seasons
fieldPoles = c(fieldFabGranites2018$pole, fieldFabGranites2016$pole)

# compute and plot shaded kamb contours (no data points will plot) for foliation
fieldKambs = plotKambContours(points = (fieldPoles), numNonAdapt = 5) + new_scale_color() + new_scale_fill()
fieldKambs

# plot data points on top of kamb contour plot
fieldFoliationEAN = plotEqualAreaNet(plot = fieldKambs, points =  list(fieldPoles), shape = 21, pointSize = 1) + annotate("text", color = "black", x = 1.15, y = -1.3, size = 4, label = paste0("N = ", length(fieldPoles)))


# Equal area net with Kamb contours for Cleland measurements
clelandPoles = c(fieldFabGranites2018[fieldFabGranites2018$SUPERSUITE == "Cleland Supersuite",]$pole, fieldFabGranites2016[fieldFabGranites2016$SUPERSUITE == "Cleland Supersuite",]$pole)
clelandKambs = plotKambContours(points = clelandPoles, numNonAdapt = 5) + new_scale_color() + new_scale_fill()          
fieldClelandEAN = plotEqualAreaNet(plot = clelandKambs, points =  list(clelandPoles), shape = 21, pointSize = 1) + annotate("text", color = "black", x = 1.15, y = -1.3, size = 4, label = paste0("N = ", length(clelandPoles)))
fieldClelandEAN


# Equal area net with Kamb contours for Emu Pool
EmuPoles = c(fieldFabGranites2018[fieldFabGranites2018$SUPERSUITE == "Emu Pool Supersuite",]$pole, fieldFabGranites2016[fieldFabGranites2016$SUPERSUITE == "Emu Pool Supersuite",]$pole)
EmuKambs = plotKambContours(points = EmuPoles, numNonAdapt = 5) + new_scale_color() + new_scale_fill()          
fieldEmuEAN = plotEqualAreaNet(plot =EmuKambs, points =  list(EmuPoles), shape = 21, pointSize = 1) + annotate("text", color = "black", x = 1.15, y = -1.3, size = 4, label = paste0("N = ", length(EmuPoles)))
fieldEmuEAN



# Equal area net with Kamb contours for Callina/Tambina 
TambinaPoles = c(fieldFabGranites2018[fieldFabGranites2018$SUPERSUITE == "Tambina Supersuite",]$pole, fieldFabGranites2016[fieldFabGranites2016$SUPERSUITE == "Tambina Supersuite",]$pole)
TambinaKambs = plotKambContours(points = TambinaPoles, numNonAdapt = 5) + new_scale_color() + new_scale_fill()          
fieldTambinaEAN = plotEqualAreaNet(plot =TambinaKambs, points =  list(TambinaPoles), shape = 21, pointSize = 1) + annotate("text", color = "black", x = 1.15, y = -1.3, size = 4, label = paste0("N = ", length(TambinaPoles)))
fieldTambinaEAN 

# Save equal area nets
ggsave("plotExports/fieldFoliationEAN.pdf", plot = fieldFoliationEAN, width = 2, height = 2, unit = "in")
ggsave("plotExports/fieldClelandEAN.pdf", plot = fieldClelandEAN, width = 2, height = 2, unit = "in")
ggsave("plotExports/fieldEmuEAN.pdf", plot = fieldEmuEAN, width = 2, height = 2, unit = "in")
ggsave("plotExports/fieldTambinaEAN.pdf", plot = fieldTambinaEAN, width = 2, height = 2, unit = "in")


dataAMS = amsDataMeans

# 
# Equal area nets and equal volume plots
plotEqualAreaNet(points = list(lapply(dataAMS$rotation, function(rot) rot[1,]),
                               lapply(dataAMS$rotation, function(rot) rot[2,]),
                               lapply(dataAMS$rotation, function(rot) rot[3,])),
                 color = list("blue", "green", "red"),
                 pointSize = 1.5,
                 shape = 21)

plotEqualVolumeEllipsoidRotations(dataAMS$rotation, pointSize = 1)


# Equal area net of magnetic lineations. All cores (grey) and station mean (black)
magLins = plotEqualAreaNet(points = list(lapply(amsData$rotation, function(s) s[1,]),
                                         lapply(amsDataMeans$rotation, function(s) s[1,])),
                           color = list("grey90", "black"),
                           shape = 15,
                           pointSize = 1.5,
                           edgeColor = list("grey90","black")) + scale_fill_manual(values = superSuiteColorsBW) + annotate("text", color = "grey40", x = -1.15, y = -1.3, size = 4, label = paste0("n = ", length(amsData$rotation))) + annotate("text", color = "black", x = -1.15, y = 1.3, size = 4, label = paste0("N = ", length(amsDataMeans$rotation)))

# Equal area net of magnetic K2. All cores (grey) and station mean (black)
magK2 = plotEqualAreaNet(points = list(lapply(amsData$rotation, function(s) s[2,]),
                                       lapply(amsDataMeans$rotation, function(s) s[2,])),
                         color = list("grey", "black"),
                         pointSize = 1.5,
                         shape = 17,
                         edgeColor = list("grey90","black")) + scale_fill_manual(values = superSuiteColorsBW) + annotate("text", color = "grey40", x = -1.15, y = -1.3, size = 4, label = paste0("n = ", length(amsData$rotation))) + annotate("text", color = "black", x = -1.15, y = 1.3, size = 4, label = paste0("N = ", length(amsDataMeans$rotation)))
magK2

# Equal area net of magnetic foliations. All cores (grey) and station mean (black)
magFols = plotEqualAreaNet(points = list(lapply(amsData$rotation, function(s) s[3,]),
                                         lapply(amsDataMeans$rotation, function(s) s[3,])),
                           color = list("grey", "black"),
                           pointSize = 1.5,
                           shape = 19,
                           edgeColor = list("grey90","black")) + scale_fill_manual(values = superSuiteColorsBW) + annotate("text", color = "grey40", x = -1.15, y = -1.3, size = 4, label = paste0("n = ", length(amsData$rotation))) + annotate("text", color = "black", x = -1.15, y = 1.3, size = 4, label = paste0("N = ", length(amsDataMeans$rotation)))
magFols


# Plot mean magnetic lineation with kamb contours
lineationKambs = plotKambContours(points = lapply(amsDataMeans$rotation,function(s) s[1,]), numNonAdapt = 5)
lineationKambs
magLinMeans = plotEqualAreaNet(plot = lineationKambs, points = list(lapply(amsDataMeans$rotation, function(s) s[1,])),
                               
                               color = list("black"),
                               shape = 22,
                               pointSize = 1.5,
                               edgeColor = list("white")) +
    annotate("text", color = "black", x = 1.15, y = -1.3, size = 4, label = paste0("N = ", length(amsDataMeans$rotation)))
magLinMeans


# Plot mean magnetic foliation with kamb contours
foliationKambs = plotKambContours(points = lapply(amsDataMeans$rotation,function(s) s[3,]), numNonAdapt = 5)
foliationKambs
magFolMeans = plotEqualAreaNet(plot = foliationKambs, points = list(lapply(amsDataMeans$rotation, function(s) s[3,])),
                               
                               color = list("black"),
                               shape = 21,
                               pointSize = 1.5,
                               edgeColor = list("white")) +
    annotate("text", color = "black", x = 1.15, y = -1.3, size = 4, label = paste0("N = ", length(amsDataMeans$rotation)))
magFolMeans


# Save plots
ggsave("plotExports/magLinsEAN.PDF", plot = magLins, height = 2.5, width = 2.5, unit = "in", useDingbats = FALSE)
ggsave("plotExports/magFolsEAN.PDF", plot = magFols, height = 2.5, width = 2.5, unit = "in", useDingbats = FALSE)
ggsave("plotExports/magK2EAN.PDF", plot = magK2, height = 2.5, width = 2.5, unit = "in", useDingbats = FALSE)
ggsave("plotExports/magLinMeansEAN.PDF", plot = magLinMeans, height = 2.5, width = 2.5, unit = "in", useDingbats = FALSE)
ggsave("plotExports/magFolMeanssEAN.PDF", plot = magFolMeans, height = 2.5, width = 2.5, unit = "in", useDingbats = FALSE)


#PLOT AMS EANs by supersuite
clelandAMS = amsDataMeans[amsDataMeans$SUPERSUITE == "Cleland Supersuite",]
clelandAMS = clelandAMS[is.na(clelandAMS$SUPERSUITE) == FALSE,]

clelandFolKambs = plotKambContours(points = lapply(clelandAMS$rotation,function(s) s[3,]), numNonAdapt = 5)
clelandLinKambs = plotKambContours(points = lapply(clelandAMS$rotation,function(s) s[1,]), numNonAdapt = 5)

clelandFolKambs
clelandLinKambs

clelandMagLins = plotEqualAreaNet(plot = clelandLinKambs, points = list(lapply(clelandAMS$rotation, function(s) s[1,])),
                                  color = list("black"),
                                  shape = 22,
                                  pointSize = 1.5,
                                  edgeColor = list("white")) + 
    annotate("text", color = "black", x = 1.15, y = -1.3, size = 4, label = paste0("N = ", length(clelandAMS$rotation)))


clelandMagLins


clelandMagFols = plotEqualAreaNet(plot = clelandFolKambs, points = list(lapply(clelandAMS$rotation, function(s) s[3,])),
                                  
                                  color = list("black"),
                                  shape = 21,
                                  pointSize = 1.5,
                                  edgeColor = list("white")) +
    annotate("text", color = "black", x = 1.15, y = -1.3, size = 4, label = paste0("N = ", length(clelandAMS$rotation)))


clelandMagFols

# Emu Pool Supersuiite
EmuAMS = amsDataMeans[amsDataMeans$SUPERSUITE == "Emu Pool Supersuite",]
EmuAMS = EmuAMS[is.na(EmuAMS$SUPERSUITE) == FALSE,]

EmuFolKambs = plotKambContours(points = lapply(EmuAMS$rotation,function(s) s[3,]), numNonAdapt = 5)
EmuLinKambs = plotKambContours(points = lapply(EmuAMS$rotation,function(s) s[1,]), numNonAdapt = 5)

EmuFolKambs
EmuLinKambs

EmuMagLins = plotEqualAreaNet(plot = EmuLinKambs, points = list(lapply(EmuAMS$rotation, function(s) s[1,])),
                              color = list("black"),
                              shape = 22,
                              pointSize = 1.5,
                              edgeColor = list("white")) + 
    annotate("text", color = "black", x = 1.15, y = -1.3, size = 4, label = paste0("N = ", length(EmuAMS$rotation)))


EmuMagLins


EmuMagFols = plotEqualAreaNet(plot = EmuFolKambs, points = list(lapply(EmuAMS$rotation, function(s) s[3,])),
                              
                              color = list("black"),
                              shape = 21,
                              pointSize = 1.5,
                              edgeColor = list("white")) +
    annotate("text", color = "black", x = 1.15, y = -1.3, size = 4, label = paste0("N = ", length(EmuAMS$rotation)))


EmuMagFols



# Tambina Supersuite
# Emu Pool Supersutie
TambinaAMS = amsDataMeans[amsDataMeans$SUPERSUITE == "Tambina Supersuite",]
TambinaAMS = TambinaAMS[is.na(TambinaAMS$SUPERSUITE) == FALSE,]

TambinaFolKambs = plotKambContours(points = lapply(TambinaAMS$rotation,function(s) s[3,]), numNonAdapt = 5)
TambinaLinKambs = plotKambContours(points = lapply(TambinaAMS$rotation,function(s) s[1,]), numNonAdapt = 5)

TambinaFolKambs
TambinaLinKambs

TambinaMagLins = plotEqualAreaNet(plot = TambinaLinKambs, points = list(lapply(TambinaAMS$rotation, function(s) s[1,])),
                                  color = list("black"),
                                  shape = 22,
                                  pointSize = 1.5,
                                  edgeColor = list("white")) + 
    annotate("text", color = "black", x = 1.15, y = -1.3, size = 4, label = paste0("N = ", length(TambinaAMS$rotation)))


TambinaMagLins


TambinaMagFols = plotEqualAreaNet(plot = TambinaFolKambs, points = list(lapply(TambinaAMS$rotation, function(s) s[3,])),
                                  
                                  color = list("black"),
                                  shape = 21,
                                  pointSize = 1.5,
                                  edgeColor = list("white")) +
    annotate("text", color = "black", x = 1.15, y = -1.3, size = 4, label = paste0("N = ", length(TambinaAMS$rotation)))


TambinaMagFols





ggsave("plotExports/clelandMagLins.pdf", plot = clelandMagLins, height = 2.5, width = 2.5, unit= "in")
ggsave("plotExports/EmuMagLins.pdf", plot = EmuMagLins, height = 2.5, width = 2.5, unit= "in")
ggsave("plotExports/TambinaMagLins.pdf", plot = TambinaMagLins, height = 2.5, width = 2.5, unit= "in")

ggsave("plotExports/clelandMagFols.pdf", plot = clelandMagFols, height = 2.5, width = 2.5, unit= "in")
ggsave("plotExports/EmuMagFols.pdf", plot = EmuMagFols, height = 2.5, width = 2.5, unit= "in")
ggsave("plotExports/TambinaMagFols.pdf", plot = TambinaMagFols, height = 2.5, width = 2.5, unit= "in")



# -------------Pj-T PLOTS ------------------------------

PjTplotAll = ggplot(amsData) + 
     # stat_density2d(bins = 10,aes(x = Pj, y = T, fill = ..level..), contour = TRUE, n = 100, geom = "polygon") +
     # scale_fill_gradient(low = "white", high = "grey2") +
     geom_point(aes(x = P, y = T), shape = 19,size = 0.25, alpha = 0.5) +
     coord_cartesian(expand = FALSE) +
     geom_hline(yintercept = 0) +
     xlim(c(1,2)) +
     ylim(c(-1,1)) +
     theme_classic()


PjTplotMeans = ggplot(amsDataMeans) + 
     # stat_density2d(bins = 10,aes(x = sapply(logA,ellJelinekP), y = sapply(logA, ellLodeNu), fill = ..level..), contour = TRUE, n = 100, geom = "polygon") +
     # scale_fill_gradient(low = "white", high = "grey2") +
     geom_point(aes(x = sapply(logA, ellJelinekP), y = sapply(logA, ellLodeNu)), size = 0.25, alpha = 0.5) +
     coord_cartesian(expand = FALSE) +
     xlab("Pj") +
     ylab("T") +
     geom_hline(yintercept = 0) +
     xlim(c(1,2)) +
     ylim(c(-1,1)) +
     theme_classic()

ggsave("plotExports/PjTplotAll.pdf", plot = PjTplotAll, height = 4, width = 6, unit = "in", useDingbats = FALSE)
ggsave("plotExports/PjTplotMeans.pdf", plot = PjTplotMeans, height = 4, width = 6, unit = "in", useDingbats = FALSE)



# bar graph of Cleland Supersuite microstructural classification

clelandMicro = microStructures[microStructures$SUPERSUITE == "Cleland Supersuite",]
clelandMicro = clelandMicro[is.na(clelandMicro$Classification) == FALSE,]

clelandMicroBar = ggplot(data = clelandMicro, aes(x = Classification)) +
    geom_bar(fill = superSuiteColorsBW[1]) + theme_classic() + coord_cartesian(expand = FALSE) + expand_limits(x = 0.5)+ theme(axis.title = element_blank(), axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + scale_x_discrete(drop = FALSE)

# bar graph of Emu Pool Supersuite microstructural classification

EmuMicro = microStructures[microStructures$SUPERSUITE == "Emu Pool Supersuite",]
EmuMicro = EmuMicro[is.na(EmuMicro$Classification) == FALSE,]

EmuMicroBar = ggplot(data = EmuMicro, aes(x = Classification)) +
    geom_bar(fill = superSuiteColorsBW[2]) + theme_classic() + coord_cartesian(expand = FALSE) + expand_limits(x = 0.5)+ theme(axis.title = element_blank(), axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + scale_x_discrete(drop = FALSE)

# bar graph of Tambina Supersuite microstructural classification

TambinaMicro = microStructures[microStructures$SUPERSUITE == "Tambina Supersuite",]
TambinaMicro = TambinaMicro[is.na(TambinaMicro$Classification) == FALSE,]

TambinaMicroBar = ggplot(data = TambinaMicro, aes(x = Classification)) +
    geom_bar(fill = superSuiteColorsBW[3]) + theme_classic() + coord_cartesian(expand = FALSE) + expand_limits(x = 0.5)+ theme(axis.title = element_blank(), axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + scale_x_discrete(drop = FALSE)

ggsave("plotExports/ClelandMicroBar.pdf", plot = clelandMicroBar, height = 1.5, width = 1.2, unit ="in")
ggsave("plotExports/EmuMicroBar.pdf", plot = EmuMicroBar, height = 1.5, width = 1.2, unit ="in")
ggsave("plotExports/TambinaMicroBar.pdf", plot = TambinaMicroBar, height = 1.5, width = 1.2, unit ="in")


Micro = microStructures
Micro = Micro[is.na(Micro$Classification) == FALSE,]

MicroBar = ggplot(data = Micro, aes(x = Classification)) +
    geom_bar(fill = superSuiteColorsBW[3]) + theme_classic() + coord_cartesian(expand = FALSE) + expand_limits(x = 0.5)+ theme(axis.title = element_blank(), axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + scale_x_discrete(drop = FALSE)
MicroBar

# bar graph of Tambina Supersuite microstructural classification
TambinaMicro = microStructures[microStructures$SUPERSUITE == "Tambina Supersuite",]
TambinaMicro = TambinaMicro[is.na(TambinaMicro$Classification) == FALSE,]

TambinaMicroBar = ggplot(data = TambinaMicro, aes(x = Classification)) +
    geom_bar(fill = superSuiteColorsBW[3]) + theme_classic() + coord_cartesian(expand = FALSE) + expand_limits(x = 0.5)+ theme(axis.title = element_blank(), axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + scale_x_discrete(drop = FALSE)

ggsave("plotExports/ClelandMicroBar.pdf", plot = clelandMicroBar, height = 1.5, width = 1.2, unit ="in")
ggsave("plotExports/EmuMicroBar.pdf", plot = EmuMicroBar, height = 1.5, width = 1.2, unit ="in")
ggsave("plotExports/TambinaMicroBar.pdf", plot = TambinaMicroBar, height = 1.5, width = 1.2, unit ="in")

