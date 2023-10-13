# Fig. 1: Maps

# TODO
# Overview map with rainforest (cf White) & refugia
# detailed map of Sangha/Likwala region (bbox PKM distribution) & vegetation satelite

library(ggplot2)
library(sf)
library(tidyr)

bb <- c(xmin = 14, xmax = 26, ymin = -6, ymax = 6)

# data ----
land10 <- rnaturalearth::ne_download(scale = 10, type = 'land', category = 'physical', returnclass = "sf")
coast10 <- rnaturalearth::ne_download(scale = 10, type = 'coastline', category = 'physical', returnclass = "sf")
rivers10 <- rnaturalearth::ne_download(scale = 10, type = "rivers_lake_centerlines", category = "physical", returnclass="sf")
lakes10 <- rnaturalearth::ne_download(scale = 10, type = "lakes", category = "physical", returnclass="sf")
boundary_lines_land10 <- rnaturalearth::ne_download(scale = 10, type = "boundary_lines_land", category = "cultural", returnclass="sf")

osm.rivers.lines <- geojsonsf::geojson_sf("gis/OSM_river_lines.geojson") %>% sf::st_crop(bb)
osm.rivers.poly <- geojsonsf::geojson_sf("gis/OSM_river_lakes_poly.geojson") %>%
  sf::st_make_valid() # %>% sf::st_union() %>% sf::st_crop(bb)
osm.coast.line <- geojsonsf::geojson_sf("gis/OSM_coast_lines.geojson") %>% sf::st_crop(bb)

sites <- data.table::fread(
  "https://raw.githubusercontent.com/dirkseidensticker/aSCAC/master/sites.csv", 
  encoding = "UTF-8") %>%
  sf::st_as_sf(coords = c("LONG", "LAT"), 
               remove = F, 
               crs = 4326, 
               na.fail = F)

pottery <- data.table::fread(
  "https://raw.githubusercontent.com/dirkseidensticker/aSCAC/master/potterygroups.csv", 
  encoding = "UTF-8") %>%
  dplyr::select(-ID, -DESCRIPTION)

c14 <- rbind(
  data.table::fread(
    "https://raw.githubusercontent.com/dirkseidensticker/aDRAC/master/aDRAC.csv", 
    encoding = "UTF-8"),
  data.table::fread(
    "https://raw.githubusercontent.com/dirkseidensticker/PikundaMunda_BatalimoMaluba_AAR/main/data/aDRAC_new.csv", 
    dec = ",", 
    encoding = "UTF-8")
)
# landcover data ----

rainforest <- geojsonsf::geojson_sf("gis/white1983.geojson") %>%
  st_set_crs(4326) %>%
  dplyr::filter(DESCRIPTIO %in% c("Anthropic landscapes",
                                  "Dry forest and thicket",
                                  "Swamp forest and mangrove",
                                  "Tropical lowland rainforest"))

refugia <- geojsonsf::geojson_sf("gis/Bremond_etal2017Fig1.geojson")

library("raster")
library("rgdal")

# landcover data from "Global Land Cover 2000 Project (GLC 2000)" https://ec.europa.eu/jrc/en/scientific-tool/global-land-cover
temp <- tempfile()
download.file("https://forobs.jrc.ec.europa.eu/data/products/glc2000/Africa_v5_Grid.zip",temp)
unzip(temp)
unlink(temp)
dpath <- "Grid/africa_v5/hdr.adf"
x <- new("GDALReadOnlyDataset", dpath)
getDriver(x)
getDriverLongName(getDriver(x))
hdr <- asSGDF_GROD(x)
hdr <- raster(hdr)
# extent: xmin,xmax,ymin,ymax
e  <- extent(5, 33, -17, 13) 
rfs <- crop(hdr, e) 

#rfs <- aggregate(rfs, fact = 2, fun = min) # reduce resolution

rfs.p <- rasterToPoints(rfs)
# Make the points a dataframe for ggplot & subset rainforest bands 1-7
rfs.bd1.7 <- data.frame(rfs.p)
rfs.bd1.7 <- subset(rfs.bd1.7, band1 >= 1 & band1 <= 7)
# Swamp forest
rfs.bd5 <- data.frame(rfs.p)
rfs.bd5 <- subset(rfs.bd5, band1 == 5)
# swamp bushland and grassland
rfs.bd17 <- data.frame(rfs.p)
rfs.bd17 <- subset(rfs.bd17, band1 == 17)
# water
rfs.bd26 <- data.frame(rfs.p)
rfs.bd26 <- subset(rfs.bd26, band1 == 26)

# plot maps ----

sites.text1 <- dplyr::filter(sites %>% dplyr::distinct(SITE, LAT, LONG), SITE %in% c(
  "Ikenge"
))

plt.map1 <- ggplot() + 
  geom_sf(data = land10, fill = "#ffebbe", color = NA) + 
  geom_sf(data = sf::st_union(rainforest), fill = "#73a788", color = NA) + 
  geom_sf(data = refugia, fill = "#478966", color = NA) + 
  geom_sf(data = coast10, size = .5, color = '#44afe3') + 
  geom_sf(data = rivers10, size = 1, color = '#44afe3') + 
  geom_sf(data = lakes10, fill = '#44afe3', color = NA) + 
  geom_sf(data = boundary_lines_land10, color = 'black', linetype = "dashed") + 
  geom_rect(aes(xmin = 16, xmax = 17.75, ymin = -1.2, ymax = 2), fill = NA, color = "red") +  
  
  #geom_point(data = c14 %>% dplyr::filter(POTTERY != '' & POTTERY != '-'), 
  #           aes(x = LONG, y = LAT), 
  #           shape = 21, fill = "grey", color = "white") + 
  #geom_point(data = sites %>% dplyr::distinct(SITE, LAT, LONG), 
  #           aes(x = LONG, y = LAT), 
  #           shape = 21, fill = "grey", color = "white") + 
  
  geom_point(data = sites.text1, 
             aes(x = LONG, y = LAT), 
             shape = 21, fill = "black", color = "white") + 
  
  ggrepel::geom_label_repel(data = sites.text1, 
                   aes(x = LONG, y = LAT, label = SITE), 
                   size = 2.5, 
                   label.padding = 0.1, min.segment.length = 0, 
                   fill = "black", color = "white") + 
  
  # COUNTRY NAMES
  shadowtext::geom_shadowtext(aes(x = 12, y = 5), label = "Cameroon", fontface  = "bold", colour = "white", size = 3) + 
  shadowtext::geom_shadowtext(aes(x = 12, y = -1), label = "Gabon", fontface  = "bold", colour = "white", size = 3) + 
  shadowtext::geom_shadowtext(aes(x = 20, y = 6), label = "Central Africa Rep.", fontface  = "bold", colour = "white", size = 3) + 
  shadowtext::geom_shadowtext(aes(x = 14, y = -3.5), label = "Rep. Congo", fontface  = "bold", colour = "white", size = 3) + 
  shadowtext::geom_shadowtext(aes(x = 24, y = -3), label = "Dem. Rep. Congo", fontface  = "bold", colour = "white", size = 3) + 
  shadowtext::geom_shadowtext(aes(x = 17.5, y = -12), label = "Angola", fontface  = "bold", colour = "white", size = 3) + 
  
  # RIVER NAMES
  shadowtext::geom_shadowtext(aes(x = 22.25, y = 2.5), label = "CONGO", colour = "white", size = 2) + 
  shadowtext::geom_shadowtext(aes(x = 19.7, y = 4.5), label = "UBANGI", colour = "white", size = 2) + 
  shadowtext::geom_shadowtext(aes(x = 25, y = 3.85), label = "UELE", colour = "white", size = 2) + 
  shadowtext::geom_shadowtext(aes(x = 15.75, y = 4), label = "KADÉÏ", colour = "white", size = 2, angle = -35) + 
  shadowtext::geom_shadowtext(aes(x = 22.75, y = -1.5), label = "TSHUAPA", colour = "white", size = 2, angle = -25) + 

  coord_sf(xlim = c(6, 32), 
           ylim = c(-16, 12)) + 
  theme_bw() + 
  theme(axis.title = element_blank())


sites.text2 <- dplyr::filter(sites %>% dplyr::distinct(SITE, LAT, LONG), SITE %in% c(
  "Pikunda",
  "Munda",
  "Ngombe",
  "Mitula",
  "Mobaka",
  "Itanga"
))

plt.map2 <- ggplot() + 
  geom_sf(data = land10, fill = "#ffebbe", color = NA) + 
  geom_raster(data = rfs.bd1.7, aes(y = y, x = x), fill = '#00734d') + 
  geom_raster(data = rfs.bd5, aes(y = y, x = x), fill = '#2b916a') + 
  geom_raster(data = subset(rfs.bd17, x < 20), aes(y = y, x = x), fill = '#54eeb7') + 
  geom_raster(data = rfs.bd26, aes(y = y, x = x), fill = '#44afe3') + 
  geom_sf(data = osm.rivers.lines, size = .5, color = '#44afe3') + 
  geom_sf(data = osm.rivers.poly, size = .5, fill = '#44afe3', color = '#44afe3') + 
  geom_sf(data = lakes10, fill = '#44afe3', color = NA) + 
  geom_sf(data = sites, shape = 21, fill = "grey", color = "white") +
  geom_sf(data = sites %>% dplyr::filter(grepl("Pikunda-Munda", POTTERY)), shape = 21, fill = "black", color = "white") +
  ggrepel::geom_label_repel(data = sites.text2, 
                   aes(x = LONG, y = LAT, label = SITE), 
                   size = 2.5, 
                   label.padding = 0.1, min.segment.length = 0, 
                   fill = "black", color = "white") + 
  
  shadowtext::geom_shadowtext(aes(x = 16.35, y = 1.5), label = "SANGHA", colour = "white", size = 2, angle = -40) + 
  shadowtext::geom_shadowtext(aes(x = 17.4, y = .55), label = "LIKWALA-AUX-HERBES", colour = "white", size = 2, angle = -55) + 
  shadowtext::geom_shadowtext(aes(x = 17.5, y = -1), label = "CONGO", colour = "white", size = 2, angle = 45) + 
  
  coord_sf(xlim = c(16, 17.75), 
           ylim = c(-1.2, 2)) + 
  scale_x_continuous(breaks = seq(16, 17.5, .5)) + 
  theme_bw() + 
  theme(axis.title = element_blank())

plt.map <- cowplot::plot_grid(
  plt.map1,
  plt.map2,
  rel_widths = c(1.5,1),
  labels = "AUTO")

ggsave("Fig_Map.pdf", width = 8, height = 5)
ggsave("output/Fig_Map.jpg", width = 8, height = 5)

