# macrotraces

source("scripts/head.R")

mt.f <- data.table::fread("input/TransGenTN_Macrotraces_Excerpt.csv", encoding = "UTF-8")

# no of vessel units studied:
mt.f %>% dplyr::distinct(OBJ) %>% nrow()

# stats per site & feature:
mt.f %>% 
  dplyr::left_join(obj, by = "ID") %>%
  dplyr::distinct(OBJ, SITE, FEAT) %>% 
  reshape2::dcast(SITE + FEAT ~ ., value.var = 'OBJ', fun.aggregate = length) %>%
  dplyr::rename(N = 3) %>%
  dplyr::mutate(PCT = N / sum(N) * 100)


# MCA

# build general abundance table:
mt.f.w <- mt.f %>%
  dplyr::mutate(PROPERTY = paste0(POSITION, "_", FEATURE)) %>%
  dplyr::select(OBJ, PROPERTY, ORIENTATION) %>%
  reshape2::dcast(OBJ ~ PROPERTY, value.var = "ORIENTATION") %>%
  tibble::column_to_rownames("OBJ")





# iteratively removing all rows and columns having at least two entries:
mt.f.w.filt <- mt.f.w %>%
  replace(., !is.na(.), "1") %>%
  replace(., is.na(.), "0") %>% 
  dplyr::mutate_if(is.character, as.numeric) %>% # reduce the entries to numerical
  quantAAR::itremove(minnumber = 2)

# selected samples:
rownames(mt.f.w.filt)

# selected features:
names(mt.f.w.filt)

# MCA ----
mt.f.mca.res <- mt.f.w %>% 
  tibble::rownames_to_column("SHERD") %>%
  dplyr::filter(SHERD %in% rownames(mt.f.w.filt)) %>%
  tibble::column_to_rownames("SHERD") %>%
  dplyr::select(names(mt.f.w.filt)) %>%
  FactoMineR::MCA(graph = F)

mt.f_set.mca.hcpc <- FactoMineR::HCPC(
  mt.f.mca.res,
  nb.clust = 0,
  graph = FALSE,
  iter.max = 100,
  min = 3,
  max = NULL)

FactoMineR::plot.MCA(mt.f.mca.res)
FactoMineR::plot.HCPC(mt.f_set.mca.hcpc, choice = "tree")
FactoMineR::plot.HCPC(mt.f_set.mca.hcpc, choice = "3D.map")

# plot
mt.f_set.mca.df <- mt.f.mca.res$ind$coord %>%
  data.frame() %>%
  tibble::rownames_to_column("OBJ")

mt.f_set.clust <- mt.f_set.mca.hcpc$data.clust %>%
  tibble::rownames_to_column("OBJ") %>%
  dplyr::select(OBJ, clust) %>%
  dplyr::left_join(mt.f_set.mca.df)

mt.f_set.MCA.res.var <- mt.f.mca.res$var$coord %>%
  data.frame() %>%
  tibble::rownames_to_column("var")

ddata <- ggdendro::dendro_data(mt.f_set.mca.hcpc$call$t$tree, type = "rectangle")

cowplot::plot_grid(
  
  ggplot(ggdendro::segment(ddata)) + 
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
    geom_text(data = ddata$labels %>% dplyr::left_join(mt.f_set.clust %>% dplyr::mutate(label = OBJ), by = "label"), 
              aes(x = x, y = -0.005, label = label, color = clust), hjust = 0) + 
    coord_flip() + 
    scale_x_continuous(limits = c(-1, nrow(ddata$labels)+1), expand = c(0, 0)) + 
    scale_y_reverse(limits = c(.3, -.3)) + 
    theme(legend.position = c(.1, .88), 
          axis.title.y = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank()),
  cowplot::plot_grid(
    ggplot(mt.f_set.clust, 
           aes(x = Dim.1, y = Dim.3, color = clust, label = OBJ)) + 
      geom_hline(yintercept = 0, linetype = "dashed") + 
      geom_vline(xintercept = 0, linetype = "dashed") + 
      geom_point() + 
      ggrepel::geom_text_repel() + 
      theme(legend.position = "none"),
    ggplot(mt.f_set.MCA.res.var, 
           aes(xend = Dim.1, x = 0, yend = Dim.3, y = 0)) + 
      geom_segment(arrow = arrow(length = unit(0.25, "cm"))) + 
      ggrepel::geom_text_repel(data = mt.f_set.MCA.res.var, 
                      aes(x = Dim.1, y = Dim.3, label = var)),
    ggplot(mt.f_set.clust, 
           aes(x = Dim.1, y = Dim.2, color = clust, label = OBJ)) + 
      geom_hline(yintercept = 0, linetype = "dashed") + 
      geom_vline(xintercept = 0, linetype = "dashed") + 
      geom_point() + 
      ggrepel::geom_text_repel() + 
      theme(legend.position = "none"),
    ggplot(mt.f_set.MCA.res.var, 
           aes(xend = Dim.1, x = 0, yend = Dim.2, y = 0)) + 
      geom_segment(arrow = arrow(length = unit(0.25, "cm"))) + 
      ggrepel::geom_text_repel(data = mt.f_set.MCA.res.var, 
                      aes(x = Dim.1, y = Dim.2, label = var)),
    nrow = 2),
  nrow = 1, rel_widths = c(1, 3)
)
#ggsave("Fig_Macrotraces.pdf", width = 16, height = 10, bg = "white")
ggsave("output/Fig_Macrotraces.jpg", width = 16, height = 10, bg = "white")

sink(file = "output/Fig_Macrotraces.jpg.txt")
mt.f_set.mca.hcpc$desc.var
sink(file = NULL)

mt.f_set.mca.hcpc$desc.var

