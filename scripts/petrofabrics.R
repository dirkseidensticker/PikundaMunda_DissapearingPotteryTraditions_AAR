# petrofarbics

source("scripts/head.R")

petro <- data.table::fread("input/TransGenTN_Petrography_Excerpt.csv", encoding = "UTF-8") %>%
  tibble::column_to_rownames("SHERD") %>%
  dplyr::mutate_all(dplyr::na_if,"")

# iteratively removing all rows and columns having at least two entries:
petro.filt <- petro %>%
  replace(., !is.na(.), "1") %>%
  replace(., is.na(.), "0") %>% 
  dplyr::mutate_if(is.character, as.numeric) %>% # reduce the entries to numerical
  quantAAR::itremove(minnumber = 2)


# selected samples:
rownames(petro.filt)

# selected features:
names(petro.filt)

petro.MCA.res <- petro %>%
  tibble::rownames_to_column("SHERD") %>%
  dplyr::filter(SHERD %in% rownames(petro.filt)) %>%
  tibble::column_to_rownames("SHERD") %>%
  dplyr::select(names(petro.filt)) %>%
  FactoMineR::MCA(graph = F)

cowplot::plot_grid(
  factoextra::fviz_mca_ind(petro.MCA.res, repel = T, label = "ind") + 
    coord_equal(),
  factoextra::fviz_mca_var(petro.MCA.res, repel = T) + 
    coord_equal(),
  nrow = 1
)
ggsave("output/Fig_Petrofarbics1.jpg", width = 12, height = 6)
