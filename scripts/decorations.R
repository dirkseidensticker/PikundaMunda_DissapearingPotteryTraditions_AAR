# Analysis of decorational patterns

# TOD
# load data from nwCongoDB
# MCA on decorations & groupings

# cf http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/114-mca-multiple-correspondence-analysis-in-r-essentials/

source("scripts/head.R")


library(RSQLite)
temp <- tempfile()
utils::download.file(
  "https://github.com/dirkseidensticker/nwCongo/raw/master/data/base/nwCongoDB.sqlite", 
  temp, 
  mode = "wb", 
  quiet = TRUE)
con <- DBI::dbConnect(SQLite(), dbname = temp)

# get data from nwCongoDB
d <- DBI::dbGetQuery(
  con,
  "SELECT 
  	t_Obj.objID AS objID,
	  [t_Ort].[ort_name] || ' (Fpl. ' || [t_ort].[Kat-Nr] || ')' AS Ort,
	  [t_Ort].[ort_kurz] || ' ' || [t_Komplex].[bef_nr] AS FEAT,
    t_K_Pos.posName,
    t_K_Verz.verzName AS Typ,
    t_K_Verz.verzBeschr,
    t_obj.Form_Gef,
    t_obj.Typ AS POTTERY
  FROM ((t_Ort LEFT JOIN t_Komplex ON t_Ort.ortID = t_Komplex.ortID)
    LEFT JOIN t_Obj ON t_Komplex.komplexID = t_Obj.komplexID) t_Obj INNER JOIN 't_ObjPosVerz' ON t_Obj.objID = 't_ObjPosVerz'.objID
	  INNER JOIN t_K_Pos ON 't_ObjPosVerz'.posID = t_K_Pos.posID
  	INNER JOIN t_K_Verz ON 't_ObjPosVerz'.verzID = t_K_Verz.verzID
  WHERE t_obj.Typ = 'PKM'")

d.feat <- d %>% dplyr::distinct(objID, FEAT)

head(d)

d <- d %>% dplyr::mutate(Typ = paste0("d", Typ))

d %>% dplyr::distinct(objID) %>% nrow()

# replace anything other than E3 & F3 with 'others' & 'indet' (no entries) to test dependecy between decorations and main vessel types
d.form.gef <- rbind(
  d %>%
    dplyr::distinct(objID, Form_Gef) %>%
    dplyr::filter(Form_Gef %in% c("E3", "F3")), 
  d %>%
    dplyr::distinct(objID, Form_Gef) %>%
    dplyr::filter(!(Form_Gef %in% c("E3", "F3"))) %>%
    dplyr::mutate(Form_Gef = "others")
)

d.w <- d %>% 
  dplyr::filter(!is.na(Typ)) %>%
  #dplyr::filter(objID != 3123) %>% # outlier in mca
  reshape2::dcast(objID ~ Typ, 
                  value.var = "POTTERY", fun.aggregate = length) %>%
  tibble::column_to_rownames("objID") %>%
  quantAAR::itremove(minnumber = 2) %>%
  dplyr::select(-"dNA") # %>%
  #tibble::rownames_to_column("objID") %>%
  #dplyr::mutate(objID = as.integer(objID)) %>%
  #dplyr::left_join(d.form.gef, by = "objID") %>%
  #tibble::column_to_rownames("objID")



d.ca.res <- d.w %>% 
  FactoMineR::CA(graph = F)


fviz_ca_biplot(d.ca.res, repel = F)

# way too much going on here; only test vessels of type E3 & F3 against each other!

d.w <- d %>% 
  dplyr::filter(!is.na(Typ) & Form_Gef %in% c("E3", "F3")) %>%
  reshape2::dcast(objID ~ Typ, 
                  value.var = "POTTERY", fun.aggregate = length) %>%
  tibble::column_to_rownames("objID") %>%
  quantAAR::itremove(minnumber = 2)
  

d.ca.res <- d.w %>% 
  FactoMineR::CA(graph = F)


fviz_ca_biplot(d.ca.res, repel = F)


d.ca.res.obj <- d.ca.res[["row"]][["coord"]] %>% 
  as.data.frame() %>%
  tibble::rownames_to_column("objID") %>%
  dplyr::mutate(objID = as.integer(objID)) %>%
  dplyr::left_join(d.form.gef, by = "objID") %>%
  dplyr::left_join(d.feat, by = "objID")
  
d.ca.res.dec <- d.ca.res[["col"]][["coord"]] %>% 
  as.data.frame()

ggplot() + 
  geom_segment(data = d.ca.res.dec, 
               aes(x = 0, xend = `Dim 1`, 
                   y = 0, yend = `Dim 2`), 
               arrow = arrow(type = "closed", 
                             length = unit(0.02, "npc")), 
               color = "grey") + 
  ggrepel::geom_text_repel(data = d.ca.res.dec %>% tibble::rownames_to_column("dec"), 
                           aes(x = `Dim 1`, y = `Dim 2`, label = dec)) + 
  geom_point(data = d.ca.res.obj, 
             aes(x = `Dim 1`, y = `Dim 2`, color = FEAT, shape = Form_Gef), size = 2) + 
  stat_ellipse(data = d.ca.res.obj, 
               aes(x = `Dim 1`, y = `Dim 2`, linetype = Form_Gef, group = Form_Gef)) + 
  coord_equal() + 
  theme_light()
ggsave("output/Fig_DecorationCA_GefForm.jpg", width = 8, height = 8)
