# spd etc main fig

# main chronology figure ----

library(concaveman)
library(ggplot2)
library(cowplot)
library(spdep)
library(sf)
library(geojsonsf)
library(tidyr)

library(rcarbon)
library(parallel)
ncores = (detectCores() - 1)

c14 <- data.table::fread(
  "https://raw.githubusercontent.com/dirkseidensticker/aDRAC/master/aDRAC.csv", 
  encoding = "UTF-8"
)

sites <- data.table::fread(
  "https://raw.githubusercontent.com/dirkseidensticker/aSCAC/master/sites.csv", 
  encoding = "UTF-8")  %>%
  sf::st_as_sf(crs = 4326, 
               coords = c("LONG", 
                          "LAT"), 
               remove = FALSE,
               na.fail = F)

pottery <- data.table::fread(
  "https://raw.githubusercontent.com/dirkseidensticker/aSCAC/master/potterygroups.csv", 
  encoding = "UTF-8")

sites.meta <- sites %>% 
  dplyr::left_join(pottery, by = "POTTERY")

regions <- geojsonsf::geojson_sf("gis/regions.geojson")

# functions ----

# spd & model testing ----

min.14C.age <- 0 # rcarbonTest fails on 0 bp dates
n.bins <- 100 # time window for binning dates (in years) within the same site	
n.simulations <- 1000 # nr of simulations
timeRange = c(3000, 0) # set the timerange of analysis in calBP, older date first
breaks = seq(3000, 0, -200) # 200 year blocks

xlim <- c(-1000,1900)

cal.norm <- FALSE
calCurve <- 'intcal20'

rcarbonModelTestValues <- function(dates, 
                                   nullhypothesis){
  
  # 3.1 calibrate & binning
  cal = rcarbon::calibrate(x = dates$C14AGE, 			
                           errors = dates$C14STD,
                           calCurves = calCurve, 
                           ncores = ncores, 
                           normalised = cal.norm)
  
  bins <- binPrep(sites = dates$SITE,
                  ages = dates$C14AGE,
                  h = n.bins # xx years cut off value
  ) 
  
  print(paste("based on", nrow(dates), "dates and", length(unique(bins)), "bins"))
  
  # 3.2 run modelTest
  print(paste(nullhypothesis, "nullhypothesis modelTest"))
  
  if(nullhypothesis == "logistic"){
    # 3.2.1 logistic model
    
    # logistic growth model
    # see https://cran.r-project.org/web/packages/rcarbon/vignettes/rcarbon.html
    
    spd = rcarbon::spd(
      cal,
      timeRange = timeRange,
      bins = bins,
      runm = kernel
    )
    
    # define xmid for nls function/inflection point of the curve as max of SPD
    logFit.mid <- plyr::round_any(
      spd$grid[which.max(spd$grid$PrDens),"calBP"],
      500, 
      f = ceiling
    )
    
    logFit <- stats::nls(
      PrDens ~ SSlogis(calBP, 
                       Asym, 
                       xmid, 
                       scale
      ),
      data = spd$grid,
      control = stats::nls.control(
        maxiter = 200
      ),
      start = list(Asym = .3, # TODO: default .2 causes error with intcal20 !!!
                   xmid = logFit.mid, # start value
                   scale = -100
      )
    )
    
    logFitDens = data.frame(
      calBP = spd$grid$calBP,
      PrDens = stats::SSlogis(
        input = spd$grid$calBP,
        Asym = coefficients(logFit)[1],
        xmid = coefficients(logFit)[2],
        scal = coefficients(logFit)[3]
      )
    )
    
    model <- rcarbon::modelTest(
      cal, 
      errors = dates$C14STD, 
      bins = bins,
      nsim = n.simulations,
      timeRange = timeRange, 
      model = "custom",
      predgrid = logFitDens, 
      runm = kernel,
      raw = TRUE
    )
    
  } else { # 3.2.2 build in models
    
    model <- modelTest(
      cal, 
      errors = dates$C14STD, 
      bins = bins, 
      nsim = n.simulations, 
      ncores = ncores, 
      timeRange = timeRange, 
      model = nullhypothesis, 
      runm = kernel, 
      raw = TRUE
    )
    
  }
  
  # 3.3 extract raw results from modelTest
  
  # 3.3.1 spd & ci ----
  model.res <- data.frame(model$result)
  model.res$calBCAD <- 1950 - model.res$calBP
  model.res$model <- nullhypothesis
  
  # 3.3.2 booms&bust ----
  plot(get("model")) # NEVER TO BE REMOVED; PLOT NEEDS TO BE DRAWN!!!
  x = plot(get("model"), 
           bbty = 'n', 
           main = "model")
  
  lst <- list()
  cnt <- 1
  # 3.3.2.1 booms ----
  if(length(x$booms) != 0) {
    for(i in 1:length(x$booms)){
      FROM <- max(x$booms[[i]][[2]])
      TO <- min(x$booms[[i]][[2]])
      rcarbon <- "positive"
      bb <- data.frame(FROM = FROM,
                       TO = TO,
                       rcarbon = rcarbon)
      lst[[cnt]] <- bb
      cnt <- cnt + 1
    }
  }
  # 3.3.2.2 busts ----
  if(length(x$busts) != 0) {
    for(i in 1:length(x$busts)){
      FROM <- max(x$busts[[i]][[2]])
      TO <- min(x$busts[[i]][[2]])
      rcarbon <- "negative"
      bb <- data.frame(FROM = FROM,
                       TO = TO,
                       rcarbon = rcarbon)
      lst[[cnt]] <- bb
      cnt <- cnt + 1
    }
  }
  bb.res <- do.call(rbind, lst)
  if(length(bb.res) != 0){
    bb.res$MODEL <- nullhypothesis
    
    # 3.4 Convert to cal. BC/AD
    
    bb.res$FROM <- 1950 - bb.res$FROM
    bb.res$TO <- 1950 - bb.res$TO
    # cut excess data: 
    bb.res$FROM[bb.res$FROM <= xlim[1]] <- xlim[1]+1
    bb.res$TO[bb.res$TO >= xlim[2]] <- xlim[2]-1
  }
  
  # 3.5 Remove boom/bust phases earlier than modeled phases
  model.start <- head(model.res[model.res$lo > 0,], 1)$calBCAD
  bb.res <- bb.res[bb.res$FROM > model.start, ]
  
  # 3.6 Return list of results
  res <- list(
    ModelRes = model.res,
    BoomsBusts = bb.res
  )
  
  return(res)
}

# TODO
# implement regional analysis: e.g. comparison Western Central Africa vs. Congo Basin
regionalAnalysis <- function(c14, 
                             sites, 
                             pottery, 
                             sites.meta){
  
  
  
}

plot(regions)

# add region filter for only the Congo rainforest regions (A-H) 

c14.sel <- c14 %>% 
  sf::st_as_sf(coords = c("LONG", "LAT"), # convert to sf
               crs = 4326, 
               remove = F, 
               na.fail = F) %>% 
  sf::st_filter(regions %>% # filter only those in region E (Sangha/Likwala-aux-Herbes)
                  dplyr::filter(id %in% LETTERS[1:8]) %>% 
                  sf::st_union()) %>%
  dplyr::filter(C14AGE > min.14C.age & # age filtering
                  C14STD > 0 & 
                  CLASS %in% c("Ia","Ib","Ic", "IIa"))


kernel <- ceiling(mean(c14.sel$C14STD)/10)*10 # time window for running mean (in years)	based mean standard error		

# apply individual nullhypotheses
c14.sel.unf <- rcarbonModelTestValues(c14.sel, "uniform")
c14.sel.lin <- rcarbonModelTestValues(c14.sel, "linear")
c14.sel.exp <- rcarbonModelTestValues(c14.sel, "exponential")
c14.sel.log <- rcarbonModelTestValues(c14.sel, "logistic")

# combine boom/bust phases
c14.sel.bb <- rbind(c14.sel.unf[["BoomsBusts"]], 
                    c14.sel.lin[["BoomsBusts"]], 
                    c14.sel.exp[["BoomsBusts"]],
                    c14.sel.log[["BoomsBusts"]])


c14.sel.bb_refined <- c14.sel.bb
c14.sel.bb_refined <- c14.sel.bb_refined [order(c14.sel.bb_refined$FROM),]

# freq of sites ----

# Input
# - Siteslist: 'chronosites' list of sites and their pottery styles; each line should contain the site name (name) & the style (label)
# - Chronology data: a list of chronological units (labels from Siteslist) and From / To dating
StyleSitesProb <- function(sites, pottery){
  # n sites
  sites.n <- as.data.frame(table(sites$POTTERY))
  colnames(sites.n) <- c("Style", "Freq")
  
  # select only style that exists in study area
  pottery <- pottery[pottery$POTTERY %in% sites$POTTERY, ]
  
  # iterate over individual style:
  style.site.freq <- NULL
  for(i in 1:nrow(pottery)){
    style <- as.character(pottery[i,"POTTERY"])
    
    # normal distribution
    # lower & upper rounded
    # TODO find a way to attribute start/end dates that fall only partially into a century (e.g. end of 1850)
    lower <- plyr::round_any(pottery[i,]$FROM, 100, f = floor)
    upper <- plyr::round_any(pottery[i,]$TO, 100, f = ceiling)
    bins <- (upper - lower) / 100
    x <- seq(1, bins, 1)
    stdabw <- function(x) {n=length(x) ; sqrt(var(x) * (n-1) / n)} # descriptive standart deviation
    y <- dnorm(x, mean(x), stdabw(x))
    y <- y * 1 / sum(y) # compensate for sum(y) around 92% - stretch to 100%
    # TODO find better distribution model
    
    # no of sites
    n <- dplyr::filter(sites.n, Style == style)$Freq
    
    # age classes
    age.class <- seq(lower+50, upper-50, 100)
    
    # in case pottery covers only one century
    if(length(age.class) < 2){
      class <- age.class
      frq <- n
      group <- style
      style.frq <- data.frame(class, frq, group)
    }else{
      # iterate over probabilities
      style.frq <- NULL
      for(j in 1:length(y)){
        frq = n * y[j]
        class <- age.class[j]
        style.frq <- rbind(style.frq, data.frame(class, frq))
      }
      style.frq$group <- style
    }
    style.site.freq <- rbind(style.site.freq, style.frq)
  }
  return(style.site.freq)
}
sites.freq <- StyleSitesProb(sites, pottery)
sites.freq.sum <- sites.freq %>% 
  dplyr::group_by(class) %>%
  dplyr::summarise(sum(frq))

# metadata on pottery styles ----

breaks <- seq(-1000, 2000, 100)
class <- seq(1,length(breaks), 1)
breaks <- data.frame(breaks, class)
for(i in 1:nrow(breaks)){
  breaks[i, "labels"] <- paste0(breaks[i,"class"], ": ", breaks[i,"breaks"], "/", breaks[i+1,"breaks"])
}
# Frequency of sites per pottery group
pottery.sites.freq <- as.data.frame(stats::aggregate(SITE ~ POTTERY, 
                                                     data =sites, 
                                                     FUN = length)) 

# Area per pottery group (Convex hull)
# see https://github.com/joelgombin/concaveman
id <- dplyr::filter(pottery.sites.freq, SITE > 3)
pottery.sites.area <- sf::st_multipolygon()
pottery.sites.area <- st_sf(polygons = st_sfc(st_polygon()))
sf::st_crs(pottery.sites.area) <- (4326)
pottery.sites.area$POTTERY <- NA
for(i in 1:nrow(id)){
  sites.f <- dplyr::filter(sites, POTTERY == id[i,1])
  conc.hull <- concaveman(sites.f)
  conc.hull$POTTERY = id[i, "POTTERY"]
  pottery.sites.area <- rbind(pottery.sites.area, conc.hull)
}
pottery.sites.area$AREA <- pottery.sites.area %>%
  sf::st_make_valid() %>%
  sf::st_area()
pottery.sites.area$AREA <- as.numeric(pottery.sites.area$AREA)/1E9 # convert m2 into k(ilo) km2
pottery.sites.area <- dplyr::filter(pottery.sites.area, AREA > 0)
pottery.sites.area <- pottery.sites.area %>%
  dplyr::left_join(pottery, by = "POTTERY")
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
pottery.cent$AGE.jitter 	<- jitter(pottery.cent$AGE, 2)
# Frequency of pottery groups per 100 years
pottery.cent.freq <- as.data.frame(table(pottery.cent$AGE))
pottery.cent.freq$Var1 <- as.numeric(as.character(pottery.cent.freq$Var1))
# merge into meta tables
pottery.cent.meta <- pottery.cent %>%
  dplyr::select(-DESCRIPTION) %>%
  dplyr::left_join(pottery.sites.freq, by = "POTTERY") %>%
  dplyr::left_join(pottery.sites.area, by = "POTTERY")
sites.cent <- merge(x = sites, 							# merge sites per style with class (200-year century list)
                    y = dplyr::select(pottery.cent, -DESCRIPTION), 
                    by = "POTTERY", 
                    allow.cartesian = TRUE)

# Distance of sites pertaining to the same style
# see https://github.com/dirkseidensticker/HumActCentralAfrica_Paper/blob/main/response_eLetter_Giresse_etal.Rmd#L95-L159
index <- unique(sites$POTTERY)
res.lst <- list()
for(i in 1:length(index)){
  sel <- dplyr::filter(sites, POTTERY == index[i])
  if(nrow(sel) >= 5){
    sel.cords <- sf::st_coordinates(sel)
    sel.knn <- spdep::knearneigh(sel.cords, 
                                 k = 4, 
                                 longlat = TRUE)
    sel.dist <- spdep::nbdists(
      spdep::knn2nb(sel.knn), 
      sel.cords, 
      longlat =  TRUE
    )
    res.lst[[i]] <- data.frame(POTTERY = index[i], 
                               MEDIAN = median(unlist(sel.dist)))
  }
}
knn.res <- do.call(rbind, res.lst)
pottery.knn <- merge(
  x = pottery, 
  y = knn.res, 
  by = "POTTERY")
breaks <- seq(-1000, 2000, 100)
class <- seq(1,length(breaks), 1)
breaks <- data.frame(breaks, class)
for(i in 1:nrow(breaks)){breaks[i, "labels"] <- paste0(breaks[i,"class"], ": ", breaks[i,"breaks"], "/", breaks[i+1,"breaks"])}
pottery.res <- data.frame(matrix(ncol = ncol(pottery.knn)+1, nrow = 0))
x <- c(names(pottery.knn), "CLASS")
colnames(pottery.res) <- x
for (i in 1:length(pottery.knn$POTTERY)){
  for (j in 1:(nrow(breaks)-1)) {
    if(pottery.knn[i,"TO"] > breaks[j,"breaks"] & 
       pottery.knn[i,"FROM"] < breaks[j+1,"breaks"]){
      l <- pottery.knn[i,]
      l$CLASS <- breaks[j,"labels"]
      pottery.res <- rbind(pottery.res, as.data.frame(l))
    }
  }
}
pottery.res$AGE <- (as.numeric(sub("/.*", "", sub(".*? ", "", pottery.res$CLASS))) + as.numeric(sub(".*/", "", sub(".*? ", "", pottery.res$CLASS)))) / 2
pottery.res$AGE.jitter 	<- jitter(pottery.res$AGE, 2)

# only region E

c14.sel.e <- c14 %>% 
  sf::st_as_sf(coords = c("LONG", "LAT"), # convert to sf
               crs = 4326, 
               remove = F, 
               na.fail = F) %>% 
  sf::st_filter(regions %>% # filter only those in region E (Sangha/Likwala-aux-Herbes)
                  dplyr::filter(id == "E") %>% 
                  sf::st_union()) %>%
  dplyr::filter(C14AGE > min.14C.age & # age filtering
                  C14STD > 0 & 
                  CLASS %in% c("Ia","Ib","Ic", "IIa"))

cal.e = rcarbon::calibrate(x = c14.sel.e$C14AGE, 			
                         errors = c14.sel.e$C14STD,
                         calCurves = calCurve, 
                         ncores = ncores, 
                         normalised = cal.norm)

bins.e <- binPrep(sites = c14.sel.e$SITE,
                ages = c14.sel.e$C14AGE,
                h = n.bins # xx years cut off value
)

spd.e = rcarbon::spd(
  cal.e,
  timeRange = timeRange,
  bins = bins.e,
  runm = kernel
)



# plots ----

ymax <- as.data.frame(c14.sel.exp[["ModelRes"]]) ; ymax <- max(ymax$PrDens)			
bb.rect.ymax <- .45

spd.test.plt <- ggplot() + 
  geom_vline(xintercept = seq(-500, 1500, 500), linetype="dashed", color = "grey") + 
  
  #geom_ribbon(data = c14.sel.unf[["ModelRes"]], 
  #            aes(x = calBCAD, ymin = lo, ymax = hi), alpha = .1) + 
  #geom_ribbon(data = c14.sel.lin[["ModelRes"]], 
  #            aes(x = calBCAD, ymin = lo, ymax = hi), alpha = .1) + 
  #geom_ribbon(data = c14.sel.exp[["ModelRes"]], 
  #            aes(x = calBCAD, ymin = lo, ymax = hi), alpha = .1) + 
  #geom_ribbon(data = c14.sel.log[["ModelRes"]], 
  #            aes(x = calBCAD, ymin = lo, ymax = hi), alpha = .1) + 
  
  geom_ribbon(data = c14.sel.log[["ModelRes"]], 
              aes(x = calBCAD, ymin = lo, ymax = hi), alpha = .3) + 
  geom_line(data = c14.sel.exp[["ModelRes"]], # spd from one of the models
            aes(x = calBCAD, y = PrDens)) + 
  
  geom_line(data = spd.e[["grid"]], 
            aes(x = 1950 - calBP, y = PrDens), 
            linetype = "dashed") + 
  
  geom_rect(data = c14.sel.bb_refined,  
            aes(xmin = FROM, xmax = TO, ymin = 0, ymax = bb.rect.ymax, fill = rcarbon), alpha = .1) + 
  scale_fill_manual(values = c("#618cff", "#f8766d")) + 
  scale_x_continuous("", 
                     limits = c(-1000, 2000), 
                     breaks=seq(-1000,2000,500), 
                     expand = c(0, 0)) + 
  scale_y_continuous("Summed probability", 
                     expand = c(0, 0)) + 
  theme_classic() + 
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())




sites.freq.plt <- ggplot(sites.freq.sum, aes(x = class, y = `sum(frq)`)) + 
  geom_vline(xintercept = seq(-500, 1500, 500), linetype="dashed", color = "grey") + 
  geom_bar(stat="identity", fill = "white", color = "black", width = 75) + 
  scale_x_continuous("", 
                     limits = c(-1000, 2000), 
                     breaks=seq(-1000,2000,500), 
                     expand = c(0, 0)) + 
  scale_y_continuous("Number of \n contemporaneous sites", 
                     expand = c(0, 0)) + 
  theme_classic() + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())

freq.plt <- ggplot() + 
  geom_vline(xintercept = seq(-500, 1500, 500), linetype="dashed", color = "grey") + 
  geom_bar(data = pottery.cent.freq, 
           aes(x = Var1, 
               weight = Freq), 
           fill = "white", 
           color = "black", 
           width = 75) + 
  scale_x_continuous("", 
                     limits = c(-1000, 2000), 
                     breaks=seq(-1000,2000,500), 
                     expand = c(0, 0)) + 
  scale_y_continuous("Number of  \n pottery groups", 
                     expand = c(0, 0)) + 
  theme_classic() + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())

qty.sites.plt <- ggplot() + 
  geom_vline(xintercept = seq(-500, 1500, 500), linetype="dashed", color = "grey") + 
  geom_boxplot(data = pottery.cent.meta, 
               aes(x = AGE, 
                   y = SITE, 
                   group = AGE), 
               outlier.shape = 3, 
               width = 75) + 
  scale_x_continuous("", 
                     limits = c(-1000, 2000), 
                     breaks=seq(-1000,2000,500), 
                     expand = c(0, 0)) + 
  scale_y_sqrt("No. of sites\n per pottery gr.", 
               expand = c(0, 0)) + 
  theme_classic() + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())

area.plt <- ggplot() + 
  geom_vline(xintercept = seq(-500, 1500, 500), linetype="dashed", color = "grey") + 
  geom_boxplot(data = pottery.cent.meta, 
               aes(x = AGE, 
                   y = AREA, 
                   group = AGE), 
               outlier.shape = 3, 
               width = 75) + 
  scale_x_continuous("", 
                     limits = c(-1000, 2000), 
                     breaks=seq(-1000,2000,500), 
                     expand = c(0, 0)) + 
  scale_y_sqrt("Distr. area of \n pottery gr. (1000 km^2)", 
               expand = c(0, 0)) + 
  theme_classic() + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())

dist.plt <- ggplot(pottery.res %>% dplyr::filter(AGE < 1800), 
                   aes(x = AGE.jitter, 
                       y = MEDIAN, 
                       group = AGE)) +  
  geom_vline(xintercept = seq(-500, 1500, 500), linetype="dashed", color = "grey") + 
  geom_boxplot(outlier.shape = 3, 
               width = 75) + 
  scale_x_continuous("cal BCE/CE", 
                     limits = c(-1000, 2000), 
                     breaks = seq(-1000, 2000, 100), 
                     expand = c(0, 0)) +
  scale_y_continuous("Median\ndistance (km)", expand = c(0, 0)) + 
  theme_classic()

plt <- cowplot::plot_grid(
  spd.test.plt, 
  sites.freq.plt,
  freq.plt,
  qty.sites.plt, 
  area.plt, 
  dist.plt, 
  ncol = 1, 
  align = "v", axis = "lr", 
  labels = "AUTO", 
  rel_heights = c(2, 1, 1, 1, 1, 1.2))

ggsave("output/Fig_SPD.jpg", plt, width = 10, height = 10)
ggsave("Fig_SPD.pdf", plt, width = 12, height = 9.75)
