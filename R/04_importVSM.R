# import VSM data
AME16_099vsm = read.csv("data/Magnetics/NR_VSM/AME16_099_vsmR.csv")
AME16_057vsm = read.csv("data/Magnetics/NR_VSM/AME16_057_vsmR.csv")
AME18_015vsm = read.csv("data/Magnetics/NR_VSM/AME18_015_vsmR.csv")
AME18_250vsm = read.csv("data/Magnetics/NR_VSM/AME18_250_vsmR.csv")
AME18_020vsm = read.csv("data/Magnetics/NR_VSM/AME18_020_vsmR.csv")
AME18_114vsm = read.csv("data/Magnetics/NR_VSM/AME18_114_vsmR.csv")
AME18_127vsm = read.csv("data/Magnetics/NR_VSM/AME18_127_vsmR.csv")
AME18_131vsm = read.csv("data/Magnetics/NR_VSM/AME18_131_vsmR.csv")
AME18_139vsm = read.csv("data/Magnetics/NR_VSM/AME18_139_vsmR.csv")
AME18_148vsm = read.csv("data/Magnetics/NR_VSM/AME18_148_vsmR.csv")

AME16_003vsm = read.csv("data/Magnetics/NR_VSM/AME16_003_vsmR.csv")
AME16_029vsm = read.csv("data/Magnetics/NR_VSM/AME16_029_vsmR.csv")
AME16_035vsm = read.csv("data/Magnetics/NR_VSM/AME16_035_vsmR.csv")
AME16_078vsm = read.csv("data/Magnetics/NR_VSM/AME16_078_vsmR.csv")
AME16_083vsm = read.csv("data/Magnetics/NR_VSM/AME16_083_vsmR.csv")
AME16_086vsm = read.csv("data/Magnetics/NR_VSM/AME16_086_vsmR.csv")
AME16_090vsm = read.csv("data/Magnetics/NR_VSM/AME16_090_vsmR.csv")
AME16_091vsm = read.csv("data/Magnetics/NR_VSM/AME16_091_vsmR.csv")
AME16_094vsm = read.csv("data/Magnetics/NR_VSM/AME16_094_vsmR.csv")
AME16_097vsm = read.csv("data/Magnetics/NR_VSM/AME16_097_vsmR.csv")
AME16_104vsm = read.csv("data/Magnetics/NR_VSM/AME16_104_vsmR.csv")
AME16_110vsm = read.csv("data/Magnetics/NR_VSM/AME16_110_vsmR.csv")
AME16_112vsm = read.csv("data/Magnetics/NR_VSM/AME16_112_vsmR.csv")
AME16_113vsm = read.csv("data/Magnetics/NR_VSM/AME16_113_vsmR.csv")
AME16_137vsm = read.csv("data/Magnetics/NR_VSM/AME16_137_vsmR.csv")
AME16_142vsm = read.csv("data/Magnetics/NR_VSM/AME16_142_vsmR.csv")
AME16_144vsm = read.csv("data/Magnetics/NR_VSM/AME16_144_vsmR.csv")
AME16_145vsm = read.csv("data/Magnetics/NR_VSM/AME16_145_vsmR.csv")
AME16_156vsm = read.csv("data/Magnetics/NR_VSM/AME16_156_vsmR.csv")
AME16_159vsm = read.csv("data/Magnetics/NR_VSM/AME16_159_vsmR.csv")
AME16_162vsm = read.csv("data/Magnetics/NR_VSM/AME16_162_vsmR.csv")
AME16_166vsm = read.csv("data/Magnetics/NR_VSM/AME16_166_vsmR.csv")
AME16_169vsm = read.csv("data/Magnetics/NR_VSM/AME16_169_vsmR.csv")
AME16_192vsm = read.csv("data/Magnetics/NR_VSM/AME16_192_vsmR.csv")
AME16_211vsm = read.csv("data/Magnetics/NR_VSM/AME16_211_vsmR.csv")



# import thermoSusc data
AME16_099cur = read.csv("data/Magnetics/NR_magminR/AME16_099_curR.csv")
AME16_057cur = read.csv("data/Magnetics/NR_magminR/AME16_057_curR.csv")
AME18_015cur = read.csv("data/Magnetics/NR_magminR/AME18_015_curR.csv")
AME18_250cur = read.csv("data/Magnetics/NR_magminR/AME18_250_curR.csv")
AME18_020cur = read.csv("data/Magnetics/NR_magminR/AME18_020_curR.csv")
AME18_114cur = read.csv("data/Magnetics/NR_magminR/AME18_114_curR.csv")
AME18_127cur = read.csv("data/Magnetics/NR_magminR/AME18_127_curR.csv")
AME18_131cur = read.csv("data/Magnetics/NR_magminR/AME18_131_curR.csv")
AME18_272cur = read.csv("data/Magnetics/NR_magminR/AME18_272_curR.csv")
AME18_148cur = read.csv("data/Magnetics/NR_magminR/AME18_148_curR.csv")
AME18_139cur = read.csv("data/Magnetics/NR_magminR/AME18_139_curR.csv")

DayPlotData = read.csv("data/Magnetics/DayPlotData_cleaned.csv")
DayPlotData = locationPairing(DayPlotData, locData, 9)

crsDecDeg = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" 
crsUTM51S = "+proj=utm +zone=51 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

# convert amsData into an sf (simple feature) object
DayPlotData = st_as_sf(DayPlotData, coords = c("easting","northing"), crs = crsDecDeg)

# convert amsData from decimal degrees to UTM 51S meters
DayPlotData = st_transform(DayPlotData, crs  = crsUTM51S)



DayPlotDataParry = read.csv("data/Magnetics/DayPlotDataParry(1980).csv")
DayPlotDataDankers = read.csv("data/Magnetics/DayPlotDataDankersSiguiera(1981).csv")
DayPlotDataHartstra = read.csv("data/Magnetics/DayPlotDataDankersHartstra(1982).csv")



                