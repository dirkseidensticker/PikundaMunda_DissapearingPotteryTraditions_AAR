# Fig. 2: Chronology & Time sliced maps

# TODO
# cf. PKM/BTM paper

library(ggplot2)
library(ggthemes)
library(ggrepel)
library(sf)
library(tidyr)

source("scripts/myfct.R")

bb <- c(xmin = 15, xmax = 19, ymin = -2, ymax = 3)

land10 <- rnaturalearth::ne_download(scale = 10, type = 'land', category = 'physical', returnclass = "sf") %>% sf::st_crop(bb)
coast10 <- rnaturalearth::ne_download(scale = 10, type = 'coastline', category = 'physical', returnclass = "sf") %>% sf::st_crop(bb)
rivers10 <- rnaturalearth::ne_download(scale = 10, type = "rivers_lake_centerlines", category = "physical", returnclass="sf") %>% sf::st_crop(bb)
lakes10 <- rnaturalearth::ne_download(scale = 10, type = "lakes", category = "physical", returnclass="sf") %>% sf::st_make_valid() %>% sf::st_crop(bb)
boundary_lines_land10 <- rnaturalearth::ne_download(scale = 10, type = "boundary_lines_land", category = "cultural", returnclass="sf") %>% sf::st_crop(bb)

osm.rivers.lines <- geojsonsf::geojson_sf("gis/OSM_river_lines.geojson") %>% sf::st_crop(bb)

sf_use_s2(FALSE)
osm.rivers.poly <- geojsonsf::geojson_sf("gis/OSM_river_lakes_poly.geojson") %>%
  sf::st_make_valid() %>% sf::st_crop(bb)
sf_use_s2(TRUE)

regions <- geojsonsf::geojson_sf("gis/regions.geojson")

sites <- data.table::fread(
  "https://raw.githubusercontent.com/dirkseidensticker/aSCAC/master/sites.csv", 
  encoding = "UTF-8") %>%
  sf::st_as_sf(coords = c("LONG", "LAT"), 
               remove = F, 
               crs = 4326, 
               na.fail = F)

sites.sangha.likwala <- sites %>% 
  dplyr::filter(REGION == "E")

pottery <- data.table::fread(
  "https://raw.githubusercontent.com/dirkseidensticker/aSCAC/master/potterygroups.csv", 
  encoding = "UTF-8") %>%
  dplyr::select(-ID, -DESCRIPTION)

breaks <- seq(-400, 2000, 200)
class <- seq(1,length(breaks), 1)
breaks <- data.frame(breaks, class)
for(i in 1:nrow(breaks)){
  breaks[i, "labels"] <- paste0(breaks[i,"class"], ": ", breaks[i,"breaks"], "/", breaks[i+1,"breaks"])
}

# Chronology ----
c14 <- rbind(
  data.table::fread(
    "https://raw.githubusercontent.com/dirkseidensticker/aDRAC/master/aDRAC.csv", 
    encoding = "UTF-8"),
  data.table::fread(
    "https://raw.githubusercontent.com/dirkseidensticker/PikundaMunda_BatalimoMaluba_AAR/main/data/aDRAC_new.csv", 
    dec = ",", 
    encoding = "UTF-8")
)

# rcarbon ####
library(rcarbon)
library(parallel)
ncores = (detectCores() - 1)


bayes <- data.table::fread("https://raw.githubusercontent.com/dirkseidensticker/PikundaMunda_BatalimoMaluba_AAR/main/tbl/tbl_bayesphases_comparison.csv", encoding = "UTF-8") %>%
  dplyr::filter(`Pottery Group` != "Ilambi") %>%
  dplyr::select(1, 3, 5) %>%
  dplyr::rename("POTTERY" = "Pottery Group", 
                "FROM" = "Bayesian Start", 
                "TO" = "Bayesian End")

pottery.sel <- rbind(
  bayes %>% 
    dplyr::left_join(pottery %>% dplyr::select(POTTERY, REGION, COL), by = "POTTERY"),
  pottery %>%
    dplyr::filter(!POTTERY %in% c(bayes %>% dplyr::pull(POTTERY))) %>%
    dplyr::select(POTTERY,FROM,TO,REGION, COL)) %>%
  dplyr::filter(POTTERY %in% c(sites.sangha.likwala %>% dplyr::filter(POTTERY != "indet") %>% dplyr::distinct(POTTERY) %>% dplyr::pull()))

# LOOP ----
# All dates/styles ----
datalist <- list()
filterlist <- list()
a.sel.list <- list()
cnt <- 1

a <- c14 %>% 
  sf::st_as_sf(coords = c("LONG", "LAT"), # convert to sf
               crs = 4326, 
               remove = F, 
               na.fail = F) %>% 
  sf::st_filter(regions %>% # filter only those in region E (Sangha/Likwala-aux-Herbes)
                      dplyr::filter(id == "E") %>% 
                      sf::st_union()) %>%
  dplyr::filter(C14AGE > 71 & 
                C14AGE < 6000 &
                (POTTERY != '' &  POTTERY != 'indet' &  POTTERY != '(indet)' & POTTERY != '-') & 
                CLASS %in% c("Ia", "Ib", "Ic", "IIc"))

#a <- dplyr::filter(c14, 
#                   C14AGE > 71 & 
#                   C14AGE < 6000 &
#                   (POTTERY != '' &  POTTERY != 'indet' &  POTTERY != '(indet)' & POTTERY != '-') & 
#                   CLASS %in% c("Ia", "Ib", "Ic", "IIc")
#)


styles <- pottery.sel %>% dplyr::pull(POTTERY) # styles in region

for(j in 1:length(styles)){
  print(paste("[", j, "/", length(styles), "] -", styles[j]))
  
  # > FILTER DATES ---- 
  d <- dplyr::filter(c14, grepl(styles[j], c14$POTTERY)) # filter for dates related to style
  #d <- dplyr::filter(a, grepl(styles[j], a$POTTERY)) # filter for dates related to style (only in region!)
  
  if (nrow(d) != 0) {
    d <- dplyr::filter(d, !grepl(paste0("\\(" , styles[j], "\\)"), d$POTTERY)) # remove cases in parantheses
    
    d <- d %>%
      dplyr::mutate(C14AGE = as.numeric(C14AGE), 
                    C14STD = as.numeric(C14STD)) %>%
      dplyr::filter(C14AGE > 70 & !is.na(C14STD))
    
    res <- rcarbonsum(d, oxcalnorm = TRUE)
    
    res[[1]]$median <- list(res[[2]])
    res[[1]]$start <- res[[3]]
    res[[1]]$style <- styles[j]
    res[[1]]$label <- paste0(res[[1]]$style, " (", nrow(d), ")")
    #res[[1]]$region <- id.lst[i]
    
    dat <- res[[1]]
    
    datalist[[cnt]] <- dat
    a.sel.list[[cnt]] <- a
    
    cnt <- cnt + 1
  }
} # end of loop through styles

styleprob <- do.call(rbind, datalist)
#sites <- do.call(rbind, a.sel.list)

styleprob.med <- unique(styleprob[, c("style", "median")])
styleprob.med <- unnest(styleprob.med, median)

names(styleprob.med)[names(styleprob.med) == "median"] <- "TO"

styleschrono <- pottery.sel %>%
  dplyr::filter(POTTERY %in% c(pottery.sel %>% dplyr::pull(POTTERY)))

## MERGE WITH STYLECHRONO
styleprob <- merge(
  x = styleprob, 
  by.x = "style", 
  y = styleschrono[,c("POTTERY", "FROM", "TO", "COL")], 
  by.y = "POTTERY", 
  sort = FALSE, all.x = TRUE)

styleprob$rel <- TRUE

# add not merged groups (what is left?)
B <- unique(styleprob$style)
A <- unique(styleschrono$POTTERY)
missing <- A[which(!A %in% B)]

style.m <- subset(styleschrono, POTTERY %in% missing) %>%
  dplyr::select(POTTERY, FROM, TO, REGION, COL)

style.m <- style.m[,c("POTTERY", "FROM", "TO", "REGION", "COL")]

names(style.m)[names(style.m) == "POTTERY"] <- "style"
names(style.m)[names(style.m) == "REGION"] <- "region"

style.m$label <- style.m$style

style.m$rel <- FALSE

style.m$grid.calBP <- NA
style.m$grid.PrDens <- NA
style.m$median <- NA
style.m$start <- NA

style.m <- style.m[,c("style", "grid.calBP", "grid.PrDens", "median", "start", "label", "FROM", "TO", "rel", "COL")]

styleprob <- rbind(styleprob, style.m)

style.m.lab <- style.m[,c("style", "FROM", "TO", "label", "COL")]

style.m.lab$mean <- (style.m.lab$FROM + style.m$TO) / 2

style.m.lab <- style.m.lab[,c("style", "mean", "label", "COL")]
style.m.lab

style.max <- styleprob %>% 
  dplyr::group_by(label) %>% 
  dplyr::slice(which.max(grid.calBP))
names(style.max)[2] <- "max"

style.min <- styleprob %>% 
  dplyr::group_by(label) %>% 
  dplyr::slice(which.min(grid.calBP))
names(style.min)[2] <- "min"

styleprob.lab <- merge(x = style.min, y = style.max[,c("style", "max")], by = "style")
styleprob.lab$mean <- (styleprob.lab$max + styleprob.lab$min) / 2

styleprob.lab <- styleprob.lab[,c("style", "mean", "label", "COL")]

styleprob.lab1 <- rbind(styleprob.lab, style.m.lab)

colnames(styleprob.lab1)[2] <- "TO"

style.box <- unique(styleprob[c("style", "FROM", "TO", "start", "rel", "COL")])
style.box

# plot chrono ----
p.chrono <- ggplot(data = style.box, 
                   aes(x = FROM, 
                       y = reorder(style, FROM), 
                       xend = TO, 
                       yend = style)) + 
  #geom_segment(aes(linetype = rel), alpha = 0) + 
  geom_segment(aes(linetype = rel, color = COL), 
               size = 3, alpha = 0.6) + 
  scale_linetype_manual(values = c("11", "solid")) + 
  scale_color_identity() + 
  ggnewscale::new_scale_color() + 
  geom_line(data = styleprob, 
            aes(x = grid.calBP,
                y = style,
                color = grid.PrDens), 
            size = 1.5) + 
  geom_point(data = styleprob.med, 
             aes(x = TO, y = style), 
             color = "black", fill = "white", shape = 21, size = 1) +   
  geom_label(aes(label = style), 
             #angle = 90, 
             hjust = 1, nudge_x = -50, 
             size = 2, label.size = NA, 
             fontface = "bold", 
             label.padding = unit(0.1, "lines")) + 
  scale_colour_gradient(low = "white", 
                        high = "black") + 
  scale_x_continuous("cal BCE/CE", 
                     limits = c(-700,2000), 
                     breaks = c(seq(-1400,1800,200), 1950), 
                     expand = c(0,0)) + 
  #scale_y_discrete(limits = rev) + 
  theme_bw() + 
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none", 
        strip.text.y = element_text(angle = 0), 
        strip.background = element_blank())


# Frequency of sites per pottery group ----
pottery.sites.freq <- as.data.frame(
  stats::aggregate(
    SITE ~ POTTERY, 
    data = sites, # %>% dplyr::filter(POTTERY %in% c(sites.ubangi %>% dplyr::pull(POTTERY) %>% unique())), 
    FUN = length)) 

# Area per pottery group (Convex hull)
# see https://github.com/joelgombin/concaveman

id <- dplyr::filter(pottery.sites.freq, SITE > 2)
conc.hull.lst <- list()
for(i in 1:nrow(id)){
  sites.f <- dplyr::filter(sites, POTTERY == id[i,1])
  conc.hull <- concaveman::concaveman(sites.f)
  conc.hull$POTTERY = id[i, "POTTERY"]
  conc.hull.lst[[i]] <- conc.hull
  #pottery.sites.area <- rbind(pottery.sites.area, conc.hull)
}
pottery.sites.area <- do.call(rbind, conc.hull.lst) %>%
  sf::st_make_valid()

# Frequency of pottery groups per bin
pottery.cent <- data.frame(matrix(ncol = ncol(pottery)+1, nrow = 0))
x <- c(names(pottery), "CLASS")
colnames(pottery.cent) <- x

for (i in 1:length(pottery$POTTERY)){
  for (j in 1:(nrow(breaks)-1)) {
    if(pottery[i,"TO"] > breaks[j,"breaks"] & 
       pottery[i,"FROM"] < breaks[j+1,"breaks"]){
      l <- pottery[i,]
      l$CLASS <- breaks[j,"labels"]
      pottery.cent <- rbind(pottery.cent, as.data.frame(l))
    }
  }
}
pottery.cent$AGE <- (as.numeric(sub("/.*", "", sub(".*? ", "", pottery.cent$CLASS))) + as.numeric(sub(".*/", "", sub(".*? ", "", pottery.cent$CLASS)))) / 2

pottery.cent.meta <- pottery.sites.area %>%
  dplyr::left_join(pottery.cent, by = "POTTERY") %>%
  dplyr::filter(!is.na(CLASS)) %>%
  sf::st_transform(crs = 32733) %>%
  sf::st_buffer(dist = 20e3) %>%
  sf::st_transform(crs = 4326)

# set missing colours to grey
pottery.cent.meta[pottery.cent.meta$COL == '',"COL"] <- "#808080"

# set label:
lbl <- unique(pottery.cent.meta[,c("AGE", "CLASS")])
lbl$CLASS <- sub(".*? ", "", lbl$CLASS)
lbl <- setNames(lbl$CLASS, lbl$AGE)

# labels for polygons
pottery.cent.meta.labs <- pottery.cent.meta %>%
  sf::st_centroid()

pottery.cent.meta.labs <- cbind(
  pottery.cent.meta.labs, 
  do.call(rbind, 
          sf::st_geometry(pottery.cent.meta.labs)) %>% 
    as_tibble() %>% 
    setNames(c("lon","lat"))) %>%
  dplyr::filter(lon > 15.5 & lon < 25 & lat > -4.25 & lat < 5)


# plot ----
p.timeslice <- ggplot(pottery.cent.meta) +
  geom_sf(data = osm.rivers.lines, linewidth = .25, color = 'grey') + 
  geom_sf(data = osm.rivers.poly, linewidth = .25, fill = 'grey', color = 'grey') + 
  geom_sf(data = lakes10, fill = 'grey', color = NA) + 
  #geom_sf(data = rivers10 %>%
  #          dplyr::filter(name == "Ubangi" | name == "Uele") %>%
  #          sf::st_crop(c(xmin = 16, xmax = 20, ymin = -1, ymax = 6)) %>%
  #          sf::st_union(), color = "black", size = 2) + 
  geom_sf(aes(fill = COL), alpha = .5) + 
  geom_label_repel(data = pottery.cent.meta.labs %>% dplyr::filter(POTTERY %in% c(sites.sangha.likwala %>% dplyr::pull(POTTERY) %>% unique())),
                   aes(x = lon, y = lat, label = POTTERY, fill = COL), 
                   size = 2.5, 
                   label.padding = 0.1, min.segment.length = 0, 
                   color = "white", fontface = "bold", 
                   max.overlaps = Inf) + 
  scale_fill_identity() + 
  scale_x_discrete(breaks = seq(16, 22, 1)) + 
  scale_y_discrete(breaks = seq(-2, 6, 1)) + 
  facet_wrap(AGE~., 
             labeller = labeller(AGE = lbl), 
             ncol = 6) + 
  coord_sf(xlim = c(16, 17.75), 
           ylim = c(-1.2, 2)) + 
  theme_few() + 
  theme(legend.position = "none", 
        axis.title = element_blank())

p <- cowplot::plot_grid(NULL, 
                        p.chrono, 
                        p.timeslice, 
                        ncol = 1, rel_heights = c(.1, 1, 1.5),
                        labels = c("A", "", "B"))

ggsave("output/Fig_SanghaLikwala_Pottery_TimeSlices.jpg", p, width = 8, height = 8)
ggsave("Fig_SanghaLikwala_Pottery_TimeSlices.pdf", p, width = 8, height = 8)
