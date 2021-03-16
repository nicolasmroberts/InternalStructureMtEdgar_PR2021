# Load shape file of east pilbara craton
shp <- sf::st_read("shapeFiles/1_500_000_Interpreted_Bedrock_Geology__2014.shp", as_tibble = FALSE)    # load Mt 

# Load geochronology shape file from the Geological Survey of Western Australia
geochron = sf::st_read("shapeFiles/GSWA_Geochronology.shp")
geochron = st_transform(geochron, crs = "+proj=utm +zone=51 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

# Load potassium content raster image
gammaK = raster("data/UnfilteredPctK.tif")
gammaKFrame = data.frame(rasterToPoints(gammaK))

#mtEdgarBbox <- c(xmin=119.75, xmax=120.40, ymin=-21.50, ymax=-20.85)
#EPTBbox <- c(xmin = 118.17, xmax = 120.86, ymin = -22.407, ymax = -20.12)
EPTBbox = data.frame(cbind(c(118.17, 120.86), c(-22.41, -20.12)))
names(EPTBbox) = c("easting", "northing")

bboxDF = data.frame(cbind(c(119.75, 120.40), c(-21.50, -20.85)))
names(bboxDF) = c("easting", "northing")

crsDecDeg = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" 
crsUTM51S = "+proj=utm +zone=51 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

# convert amsData into an sf (simple feature) object
bboxSF = st_as_sf(bboxDF, coords = c("easting","northing"), crs = crsDecDeg)
EPTBboxSF = st_as_sf(EPTBbox, coords = c("easting","northing"), crs = crsDecDeg)

# convert amsData from decimal degrees to UTM 51S meters
bboxSF = st_transform(bboxSF, crs  = crsUTM51S)
EPTBboxSF =st_transform(EPTBboxSF, crs  = crsUTM51S)

UTMxlim = c(bboxSF$geometry[[1]][1],bboxSF$geometry[[2]][1])
UTMylim = c(bboxSF$geometry[[1]][2],bboxSF$geometry[[2]][2])


mtEdgarBboxUTM = c(xmin =bboxSF$geometry[[1]][1], xmax = bboxSF$geometry[[2]][1], ymin = bboxSF$geometry[[1]][2],ymax = bboxSF$geometry[[2]][2] )
EPTBboxUTM = c(xmin =EPTBboxSF$geometry[[1]][1], xmax = EPTBboxSF$geometry[[2]][1], ymin = EPTBboxSF$geometry[[1]][2],ymax = EPTBboxSF$geometry[[2]][2] )

pilbaraMap = st_transform(shp, crs  = "+proj=utm +zone=51 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

mtEdgarMap = st_crop(st_buffer(pilbaraMap, 0), mtEdgarBboxUTM)
eastPilbaraMap = st_crop(st_buffer(pilbaraMap, 0), EPTBboxUTM)
geochron = st_crop(geochron, mtEdgarBboxUTM)

# colors for map units
pilbaraColors = c("Apex Basalt"                 = "#d8ebd1",
                  "Bishop Creek Monzogranite"   = "#edb8d4",
                  "Black Range Dolerite Suite"  = "#000000",
                  "Bookargemoona Tonalite"      = "#c967a8",
                  "Bridget Suite"               = "#ee94be",
                  "Budjan Creek Formation"      = "#fce59b",
                  "Bullgarina Monzogranite"     = "#ecbad6",
                  "Callina Supersuite"          = "#d13d95",
                  "Campbell Well Granodiorite"  = "#e190bd",
                  "Carbana Monzogranite"        = "#e092be",
                  "Charteris Basalt"            = "#b4d9a3",
                  "Chessman Granodiorite"       = "#e190bd",
                  "Chimingadgi Trondhjemite"    = "#e093bf",
                  "Cleaverville Formation"      = "#6ec0eb",
                  "Coonieena Basalt"            = "#c9cb9f",
                  "Coppin Gap Granodiorite"     = "#e08fbc",
                  "Cotton Well Granodiorite"    = "#e090bd",
                  "Dalton Suite"                = "#7dc79b",
                  "Davitt Syenogranite"         = "#df90bd",
                  "Duffer Formation"            = "#fad453",
                  "Emu Pool Supersuite"         = "#e190bd",
                  "Euro Basalt"                 = "#b2d9a3",
                  "Farrel Quartzite"            = "#b2b271",
                  "Fig Tree Gneiss"             = "#d465a6",
                  "Fortescue Group"             = "#d29fb0",
                  "Gap Intrusion"               = "#b2aad4",
                  "Gobbos Granodiorite"         = "#e190bd",
                  "Gorge Creek Group"           = "#cacca1",
                  "Hardey Formation"            = "#d29fb0",
                  "Homeward Bound Granite"      = "#cc3f96",
                  "Jeerinah Formation"          = "#d29fb0",
                  "Jenkin Granodiorite"         = "#df92be",
                  "Johansen Monzogranite"       = "#de92be",
                  "Joorina Granodiorite"        = "#e090bd",
                  "Kennell Granodiorite"        = "#e090bd",
                  "Kylena Formation"            = "#d29fb0",
                  "Lady Adelaide Orthogneiss"   = "#ba79df",
                  "Maddina Formation"           = "#d29fb0",
                  "McPhee Formation"            = "#80c247",
                  "Mondana Monzogranite"        = "#a773ca",
                  "Moolyella Monzogranite"      = "#f16e73",
                  "Mount Ada Basalt"            = "#b4d9a1",
                  "Mount Roe Basalt"            = "#d59eb1",
                  "Mullugunya Granodiorite"     = "#e190bd",
                  "Munganbrina Monzogranite"    = "#e090bd",
                  "Nandingarra Granodiorite"    = "#e190bd",
                  "Nob Well Intrusion"          = "#3a8fde",
                  "North Star Basalt"           = "#80c243",
                  "Owens Gully Diorite"         = "#d03d95",
                  "Panorama Formation"          = "#fcd44b",
                  "Pilbara Supergroup"          = "#948ec4",
                  "Strelley Pool Formation"     = "#cdeafa",
                  "Strutton Intrusion"          = "#a6d7bd",
                  "Tambina Supersuite"          = "#d565a6",
                  "Tumbiana Formation"          = "#d29fb0",
                  "Underwood Gneiss"            = "#cf3d95",
                  "Walgunya Trondhjemite"       = "#e190bd",
                  "Warrawoona Group"            = "#c9e2ab",
                  "Wilina Granodiorite"         = "#e08fbc",
                  "Wolline Monzogranite"        = "#ecbad6",
                  "Wyman Formation"             = "#fcd177",
                  "Zulu Granodiorite"           = "#e190bd")


superSuiteColors = c("Cleland Supersuite"   = "#edb8d4",
                  "Emu Pool Supersuite"        = "#e090bd",
                  "Tambina Supersuite"          = "#d565a6",
                  "Callina Supersuite"          = "#d13d95",
                  "Split Rock Supersuite"      = "#f16e73")


superSuiteColorsBW = c("Cleland Supersuite"   = "gray75",
                     "Emu Pool Supersuite"        = "grey45",
                     "Tambina Supersuite"          = "grey15",
                     "Callina Supersuite"          = "grey15",
                     "Split Rock Supersuite"      = "grey95")

