library(tidyverse)

# sites dating to between 600 and 1000 CE in the western Congo Basin?
# cf. Fig. S4 from Seidensticker et al. 2021 (SciAdv)

c14 <- data.table::fread("https://raw.githubusercontent.com/dirkseidensticker/aDRAC/refs/heads/master/aDRAC.csv", encoding = "UTF-8") %>%
  dplyr::mutate(
    C14AGE = as.numeric(C14AGE), 
    C14STD = as.numeric(C14STD), 
    LAT = as.numeric(LAT), 
    LONG = as.numeric(LONG))

c14.wCB <- c14[which(sf::st_intersects(
  geojsonsf::geojson_sf("gis/regions.geojson") %>%
    dplyr::filter(id == "E") %>%
    sf::st_union(), 
  c14 %>%
    sf::st_as_sf(coords = c("LONG", "LAT"), crs = 4326, remove = F, na.fail = F) , 
  sparse = FALSE)), ]

# kalibieren und auf Zeitfenster 600-100 cal CE filtern:
c14.wCB.cal <- c14.wCB %>% 
  dplyr::mutate(c14age = C14AGE,
                c14std = C14STD) %>%
  c14bazAAR::as.c14_date_list() %>%
  c14bazAAR::calibrate(choices = "calrange") %>% 
  tidyr::unnest(cols = c("calrange"))

c14.wCB.hiatus <- c14.wCB.cal %>%
  dplyr::filter(from < 1950 - 600 & to > 1950 - 1000) %>%
  dplyr::select(-dens, -from, -to, -sigma) %>%
  dplyr::distinct()

c14.wCB.hiatus.cal <- c14.wCB.hiatus %>%
  c14bazAAR::as.c14_date_list() %>%
  c14bazAAR::calibrate(choices = "calprobdistr") %>%
  tidyr::unnest(cols = c("calprobdistr")) %>%
  dplyr::mutate(SOURCE = replace(SOURCE, SOURCE == "Fay 1997; Garcin et al. 2018", "Fay 1997: 221 Tab. 6.1")) %>%
  dplyr::mutate(SOURCE = replace(SOURCE, SOURCE == "Seidensticker 2017", "Seidensticker 2021: 355-356 Appendix 2")) %>%
  dplyr::mutate(LABEL = paste0(SITE, " (", SOURCE, ")"))
  
c14.wCB.hiatus.sites <- c14.wCB.hiatus.cal %>% 
  dplyr::distinct(LABEL, LAT) %>%
  dplyr::arrange(desc(LAT)) %>%
  dplyr::pull(LABEL)

c14.wCB.hiatus.cal %>%
  dplyr::mutate(LABEL = factor(LABEL, levels = c14.wCB.hiatus.sites)) %>%
  ggplot() + 
  annotate("rect", xmin = 600, xmax = 1000, ymin = -Inf, ymax = Inf, alpha = .1, fill = "red") + 
  ggridges::geom_ridgeline(
    aes(
      x = -calage + 1950, 
      y = LABNR, 
      height = density),
    scale = 100) + 
  facet_wrap(LABEL ~ ., 
             ncol = 1,
             scales = "free_y") + 
  scale_x_continuous("yrs cal CE", 
                     breaks = seq(0, 2000, 200)) + 
  ggthemes::theme_base() + 
  theme(plot.background = element_rect(color = NA), 
        axis.title.y = element_blank())
ggsave("output/rev_fig_c14_wCB_hiatus.pdf", width = 10, height = 14)

