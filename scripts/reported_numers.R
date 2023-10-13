# basic numbers listed in the text

library(sf)
library(tidyr)

sites <- read.csv(
  "https://raw.githubusercontent.com/dirkseidensticker/aSCAC/master/sites.csv",
  encoding = "UTF-8") %>%
  st_as_sf(crs = 4326, 
           coords = c("LONG", 
                      "LAT"), 
           remove = FALSE, 
           na.fail = F)

# sites with PKM pottery
sites %>% 
  dplyr::filter(POTTERY == "Pikunda-Munda") %>% 
  dplyr::distinct(SITE)

# vessels in study
obj <- readxl::read_xlsx(
  "../../0_2021-2024 Trans-Generational Training Networks (Gent)/Samples.xlsx",
  sheet = 1
)

obj %>%
  dplyr::filter(POTTERY == "Pikunda-Munda") %>%
  dplyr::distinct(SITE, FEAT, IND) %>% 
  reshape2::dcast(SITE + FEAT ~ ., value.var = 'IND', fun.aggregate = length) %>%
  dplyr::rename(N = 3) %>%
  dplyr::mutate(PCT = N / sum(N) * 100)


# thin-sections in study:
obj %>%
  dplyr::filter(POTTERY == "Pikunda-Munda" & !is.na(tsID) & tsID < 95)

obj %>%
  dplyr::filter(POTTERY != "Pikunda-Munda" & !is.na(tsID) & tsID < 95 & (SITE == "Munda" | SITE == "Pikunda"))

# vessels macro-traces:
feat <- readxl::read_xlsx(
  "../../0_2021-2024 Trans-Generational Training Networks (Gent)/Macrotraces.xlsx",
  sheet = 1)

feat %>%
  dplyr::left_join(obj %>%
                     dplyr::select(ID, SITE, FEAT, POTTERY, TYPE), 
                   by = "ID") %>%
  dplyr::filter(POTTERY == "Pikunda-Munda") %>% # filter only Pikunda-Munda style vessels
  dplyr::select(-POTTERY) %>%
  dplyr::distinct(ID) %>% nrow()


# vessels Xrayed:
obj %>%
  dplyr::filter(Xray == "x",
                POTTERY == "Pikunda-Munda") %>%
  dplyr::distinct(SITE, FEAT, IND) %>% 
  reshape2::dcast(SITE + FEAT ~ ., value.var = 'IND', fun.aggregate = length) %>%
  dplyr::rename(N = 3) %>%
  dplyr::mutate(PCT = N / sum(N) * 100)
