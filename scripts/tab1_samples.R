# Table: Overview of all studied samples

d <- data.table::fread("input/TransGenTN_Samples_Excerpt.csv") %>%
  dplyr::left_join(data.table::fread("input/TransGenTN_Macrotraces_Excerpt.csv") %>%
                     dplyr::select(ID) %>%
                     dplyr::distinct(ID) %>%
                     dplyr::mutate(MacroTraces = "x")) %>%
  dplyr::mutate(Sample = paste0(CODE, " ", FEAT, ":", IND), 
                TYPE = dplyr::recode(TYPE, 
                                     'W' = 'Wall', 
                                     'R' = 'Rim', 
                                     'B' = 'Base', 
                                     'G' = 'Vessel')) %>%
  dplyr::mutate(tsID = dplyr::case_when(tsID >= 95 ~ NA, TRUE ~ tsID))

d <- d %>%
  dplyr::filter(!(ID %in% (d %>% dplyr::filter(is.na(tsID) & is.na(Xray) & is.na(MacroTraces)) %>% dplyr::pull(ID)))) %>%
  dplyr::mutate(TYPE = replace(TYPE, is.na(TYPE), ""),
                tsID = replace(tsID, !is.na(tsID), "x"),
                Xray = replace(Xray, is.na(Xray), ""), 
                MacroTraces = replace(MacroTraces, is.na(MacroTraces), "")) %>%
  dplyr::mutate(tsID = replace(tsID, is.na(tsID), "")) %>%
  dplyr::rename(Site = SITE, 
                Style = POTTERY, 
                Sherd = TYPE,
                Petrogr = tsID) %>%
  dplyr::select(Site, Sample, Style, Sherd, MacroTraces, Xray, Petrogr) %>%
  dplyr::arrange(Sample)

xlsx::write.xlsx2(d, "output/Tab_Samples.xlsx", row.names = F)
