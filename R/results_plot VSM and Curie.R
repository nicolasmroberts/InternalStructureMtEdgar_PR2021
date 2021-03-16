plotHysteresisLoop = function(vsmDataFrame, parasub = TRUE, addTitle = FALSE) {
     p = ggplot(data = vsmDataFrame) +
          geom_point(aes(x = Field, y = Moment), size = 0.1, color = "grey") +
          geom_hline(yintercept = 0) +
          geom_vline(xintercept = 0) +
          #xlab(element_blank()) +
          #ylab(element_blank()) +
          ggtitle(deparse(substitute(vsmDataFrame))) +
          theme_classic()    
     if (parasub == TRUE) {
          p = p + geom_point(aes(x = Field, y = MomentParasub), size = 0.1, color = "black")   
     }
     if(addTitle == TRUE){
        p = p + ggtitle(deparse(substitute(vsmDataFrame)))
     }
     p
}

plotThermoSusc = function(curDataFrame, plotCurie = TRUE) {
     
        # curDataFrame$TSUSC = curDataFrame$TSUSC/curDataFrame$TSUSC[1]
        thresh = which(curDataFrame$TEMP == max(curDataFrame$TEMP))
     curHeating = curDataFrame[1:thresh,]
     curCooling = curDataFrame[thresh+1:length(curDataFrame$TEMP),]
     curPoint = curHeating[curHeating$TSUSC == max(curHeating$TSUSC),]
     cp = curPoint$TEMP
    
     
     
     p = ggplot() +
             geom_path(data = curCooling, aes(x = TEMP, y = TSUSC), size = 0.5, color = "grey75") +
             geom_point(data = curCooling, aes(x = TEMP, y = TSUSC), size = .5, color = "grey75") +
          geom_path(data = curHeating, aes(x = TEMP, y = TSUSC), size = 0.5, color = "black") +
          geom_point(data = curHeating, aes(x = TEMP, y = TSUSC), size = .5, color = "black") +
          
          xlab("Temperature (°C)") +
          ylab("Susceptibility (SI)") +
          ggtitle(deparse(substitute(curDataFrame))) +
          theme_classic()    
     if (plotCurie == TRUE) {
          p = p +   
               geom_vline(xintercept = cp) +
               geom_text(aes(cp, min(curDataFrame$TSUSC), hjust = 1.1, vjust = -1, label = paste0(cp,"°C"))) 
     }
     p
}

HysteresisPlots = list(plotHysteresisLoop(AME16_099vsm, parasub = TRUE, addTitle = TRUE),
             plotHysteresisLoop(AME16_057vsm, parasub = TRUE,addTitle = TRUE),
             plotHysteresisLoop(AME18_015vsm, parasub = FALSE,addTitle = TRUE),
             plotHysteresisLoop(AME18_250vsm, parasub = FALSE,addTitle = TRUE),
             plotHysteresisLoop(AME18_020vsm, parasub = TRUE,addTitle = TRUE),
             plotHysteresisLoop(AME18_114vsm, parasub = FALSE,addTitle = TRUE),
             plotHysteresisLoop(AME18_127vsm, parasub = FALSE,addTitle = TRUE),
             plotHysteresisLoop(AME18_131vsm, parasub = TRUE,addTitle = TRUE),
             plotHysteresisLoop(AME18_148vsm, parasub = TRUE,addTitle = TRUE),
             plotHysteresisLoop(AME18_139vsm, parasub = TRUE,addTitle = TRUE)
)


HysteresisPlotNames = c("AME16_099vsm",
                        "AME16_057vsm",
                        "AME18_015vsm",
                        "AME18_250vsm",
                        "AME18_020vsm",
                        "AME18_114vsm",
                        "AME18_127vsm",
                        "AME18_131vsm",
                        "AME18_148vsm",
                        "AME18_139vsm")

for(i in 1:length(HysteresisPlotNames)) {
    ggsave(paste0("plotExports/",HysteresisPlotNames[i],"_Hysteresis.pdf"), plot = HysteresisPlots[[i]], width = 2, height = 2, unit = "in", useDingbats = FALSE)   
}






ThermoSuscPlots = list(plotThermoSusc(AME16_099cur, plotCurie = TRUE),
                       plotThermoSusc(AME16_057cur, plotCurie = FALSE),
                       plotThermoSusc(AME18_015cur, plotCurie = FALSE),
                       plotThermoSusc(AME18_250cur, plotCurie = FALSE),
                       plotThermoSusc(AME18_020cur, plotCurie = FALSE),
                       plotThermoSusc(AME18_114cur, plotCurie = TRUE),
                       plotThermoSusc(AME18_127cur, plotCurie = FALSE),
                       plotThermoSusc(AME18_131cur, plotCurie = FALSE),
                       plotThermoSusc(AME18_148cur, plotCurie = TRUE),
                       plotThermoSusc(AME18_139cur, plotCurie = TRUE),
                       plotThermoSusc(AME18_272cur, plotCurie = FALSE)
)

ThermoSuscPlotNames = list("AME16_099cur",
                           "AME16_057cur",
                           "AME18_015cur",
                           "AME18_250cur",
                           "AME18_020cur",
                           "AME18_114cur",
                           "AME18_127cur",
                           "AME18_131cur",
                           "AME18_148cur",
                           "AME18_139cur",
                           "AME18_272cur")

for(i in 1:length(HysteresisPlotNames)) {
    ggsave(paste0("plotExports/",ThermoSuscPlotNames[i],"_ThermoMag.pdf"), plot = ThermoSuscPlots[[i]], width = 3.5, height = 2, unit = "in", useDingbats = FALSE)   
}
thermLayout = rbind(c(1,2),
                    c(3,4),
                    c(5,6))
                    

ggsave("Appendices/ThermoSuscPlots.pdf", plot = marrangeGrob(grobs = ThermoSuscPlots, layout_matrix = thermLayout,top = NULL,
                                                             heights = unit(c(3,3,3), c("in","in","in")), 
                                                             widths = unit(c(3.25,3.25), c("in","in","in"))
       ), 
       height = 11, width = 8.5, unit = "in")


curiePlots = list(plotThermoSusc(AME18_139cur, plotCurie = TRUE) + xlab(label = NULL) + ylab(label = NULL)+ xlim(c(0,650)),
                  plotThermoSusc(AME18_148cur, plotCurie = TRUE) + xlab(label = NULL) + ylab(label = NULL)+ xlim(c(0,650)),
                  plotThermoSusc(AME18_114cur, plotCurie = TRUE) + xlab(label = NULL) + ylab(label = NULL)+ xlim(c(0,650)),
                  plotThermoSusc(AME16_057cur, plotCurie = FALSE) +
                        geom_vline(xintercept = 593.1) +
                        geom_text(aes(593.1, min(-185), hjust = 1.1, vjust = -1, label = "593.1°C")) + xlab(label = NULL) + ylab(label = NULL)+ xlim(c(0,650)),
                  plotThermoSusc(AME16_099cur, plotCurie = TRUE) + xlab(label = NULL) + ylab(label = NULL)+ xlim(c(0,650)))

curiePlotNames = c("AME18-139","AME18-148","AME18-114","AME15-057","AME18-099")

for(i in 1:5) {
        ggsave(paste0("plotExports/",curiePlotNames[i],"_curiePlotFINETUNED.pdf"), plot = curiePlots[[i]], width = 3.5, height = 2, unit = "in", useDingbats = FALSE)   
}





layout = rbind(c(1,2),
               c(1,2),
               c(1,2),
               c(1,2),
               c(1,2),
               c(1,2),
               c(1,2),
               c(1,2),
               c(1,2),
               c(1,2),
               c(1,2))


ggsave("plotExports/Hysteresis_Curie_Dallas.pdf", plot = marrangeGrob(plots, layout_matrix = layout, top = NULL), height = 3, width = 6, unit = "in")


plots_Minnesota = list(plotHysteresisLoop(AME16_003vsm,addTitle = TRUE),
                       plotHysteresisLoop(AME16_029vsm,addTitle = TRUE),
                       plotHysteresisLoop(AME16_035vsm,addTitle = TRUE),
                       plotHysteresisLoop(AME16_078vsm,addTitle = TRUE),
                       plotHysteresisLoop(AME16_083vsm,addTitle = TRUE),
                       plotHysteresisLoop(AME16_086vsm,addTitle = TRUE),
                       plotHysteresisLoop(AME16_090vsm,addTitle = TRUE),
                       plotHysteresisLoop(AME16_091vsm,addTitle = TRUE),
                       plotHysteresisLoop(AME16_094vsm,addTitle = TRUE),
                       plotHysteresisLoop(AME16_097vsm,addTitle = TRUE),
                       plotHysteresisLoop(AME16_104vsm,addTitle = TRUE),
                       plotHysteresisLoop(AME16_110vsm,addTitle = TRUE),
                       plotHysteresisLoop(AME16_112vsm,addTitle = TRUE),
                       plotHysteresisLoop(AME16_113vsm,addTitle = TRUE),
                       plotHysteresisLoop(AME16_137vsm,addTitle = TRUE),
                       plotHysteresisLoop(AME16_142vsm,addTitle = TRUE),
                       plotHysteresisLoop(AME16_144vsm,addTitle = TRUE),
                       plotHysteresisLoop(AME16_145vsm,addTitle = TRUE),
                       plotHysteresisLoop(AME16_156vsm,addTitle = TRUE),
                       plotHysteresisLoop(AME16_159vsm,addTitle = TRUE),
                       plotHysteresisLoop(AME16_162vsm,addTitle = TRUE),
                       plotHysteresisLoop(AME16_166vsm,addTitle = TRUE),
                       plotHysteresisLoop(AME16_169vsm,addTitle = TRUE),
                       plotHysteresisLoop(AME16_192vsm,addTitle = TRUE),
                       plotHysteresisLoop(AME16_211vsm,addTitle = TRUE))

layoutMatrixMinnesota = rbind(c(3,2),
                              c(3,2),
                              c(3,2),
                              c(3,2),
                              c(1))


allHysteresisPlots = c(HysteresisPlots, plots_Minnesota)


HysteresisPlotNames = c(HysteresisPlotNames, 
                        "AME16_003vsm",
                        "AME16_029vsm",
                        "AME16_035vsm",
                        "AME16_078vsm",
                        "AME16_083vsm",
                        "AME16_086vsm",
                        "AME16_090vsm",
                        "AME16_091vsm",
                        "AME16_094vsm",
                        "AME16_097vsm",
                        "AME16_104vsm",
                        "AME16_110vsm",
                        "AME16_112vsm",
                        "AME16_113vsm",
                        "AME16_137vsm",
                        "AME16_142vsm",
                        "AME16_144vsm",
                        "AME16_145vsm",
                        "AME16_156vsm",
                        "AME16_159vsm",
                        "AME16_162vsm",
                        "AME16_166vsm",
                        "AME16_169vsm",
                        "AME16_192vsm",
                        "AME16_211vsm")


HysteresisPlotNames



ggsave("plotExports/allHysteresisPlots.pdf", plot = marrangeGrob(allHysteresisPlots, nrow  = 3, ncol = 2,top = NULL), height = 9, width = 6.5, unit = "in")


thermLayout = rbind(c(1,2),
                    c(3,4),
                    c(5,6))

ggsave("plotExports/allHysteresisPlots.pdf", plot = marrangeGrob(grobs = allHysteresisPlots, layout_matrix = thermLayout,top = NULL,
                                                             heights = unit(c(3,3,3), c("in","in","in")), 
                                                             widths = unit(c(3.25,3.25), c("in","in","in"))
), 
height = 11, width = 8.5, unit = "in")

ggsave("plotExports/Hysteresis_Curie_Dallas.pdf", plot = marrangeGrob(plots, layout_matrix = layout, top = NULL), height = 3, width = 6, unit = "in")

ggsave("plotExports/Hysteresis_Minnesota.pdf", plot = marrangeGrob(plots_Minnesota, nrow = 3, ncol = 2, top = NULL), height = 11, width = 8.5 , unit = "in")


DayPlotData$SUPERSUITE = factor(DayPlotData$SUPERSUITE, levels = c("Cleland Supersuite", "Emu Pool Supersuite", "Tambina Supersuite","Callina Supersuite"))


dayPlot = ggplot(DayPlotData) +
     geom_hline(yintercept = 0.5) +
     geom_hline(yintercept = 0.02) +
     geom_vline(xintercept = 2) +
     geom_vline(xintercept = 5) +
     #coord_trans(x = "log10", y = "log10") +
        scale_x_log10(limits = c(1,10^1.5), breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x))) +
        scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x))) +
        annotation_logticks() +
        

        #geom_point(data = DayPlotDataParry, aes(x = HrHc, y = MrsMs), shape = 0) +
        geom_point(data = DayPlotDataDankers, aes(x = HrHc, y = MrsMs), shape = 1, color = "grey") +
        geom_point(data= DayPlotDataHartstra, aes(x=HcrHc, y = MrsMs), shape = 3, color = "grey") +
        xlab("Coercivity Ratio Hcr / Hc") +
        ylab("Remanence Ratio Mrs / Ms") +
        geom_point(aes(x=Hcr/Hc, y = Mrs/Ms, fill = SUPERSUITE), shape = 21, size = 2) +
        scale_fill_manual(values = superSuiteColorsBW) +
     theme_classic() +
    theme(panel.grid.major = element_blank(), 
          legend.position = c(1,1),
          legend.margin = margin(0.005,.005,.005,.005, "npc"),
          legend.justification = c(1,1),
          legend.key.height = unit(0.04, "npc"), 
          legend.direction = "vertical",
          legend.background = element_rect(fill = "white", color = "black", size = .5),
          legend.title = element_blank(),
          panel.border = element_rect(size = 1, color = "black", fill = NA))

dayPlot

ggsave("plotExports/DayPlot.pdf", plot = dayPlot, height = 4, width = 3, unit = "in", useDingbats = FALSE)
