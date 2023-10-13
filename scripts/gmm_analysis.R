# GMM Analysis

# TODO
# Eliptical fourrier transformation & harmonics
# PCA

library(ggplot2)
library(Momocs)

obj <- data.table::fread(
  "https://raw.githubusercontent.com/dirkseidensticker/aSCAC/master/potterydrawings.csv",
  encoding = "UTF-8", na.strings=c("")) %>%
  dplyr::mutate(OBJ = SOURCE)
head(obj)


# GMM ----
files <- paste0("output/gmm_masks//", 
                list.files("output/gmm_masks/"))

x <- Momocs::import_jpg(files)

ves <- Out(x) # building out object

# plot outlines:
jpeg(filename = "output/Fig_GMM_panel.jpg", width = 8, height = 8, units = "in", res = 300)
panel(ves)
dev.off()

jpeg(filename = "output/Fig_GMM_stack.jpg", width = 8, height = 8, units = "in", res = 300)
ves %>% coo_center %>% stack()
dev.off()

# Morphometrics ----

# Outline Analysis with Elliptical Fourier Transforms
coo_oscillo(ves[1], "efourier")

# Fitting Ptolemaic ellipses on the plane
Ptolemy(ves[1])

calibrate_harmonicpower_efourier(ves)
calibrate_deviations_efourier(ves)
calibrate_reconstructions_efourier(ves)

ves.f <- efourier(ves)
ves.f

#hist(ves.f, drop=0)
boxplot(ves.f, drop=1)

ves.p <- PCA(ves.f)
class(ves.p) # a PCA object, let's plot it
plot(ves.p)

# export the PCA eigenvectors
gmm_potplot_pca <- ves.p[["x"]] %>%
  data.frame() %>%
  tibble::rownames_to_column("OBJ")

write.csv(gmm_potplot_pca, "output/gmm_potplot_pca.csv", row.names = F)

library("ggimage")

gmm_potplot_pca %>% 
  dplyr::mutate(label = gsub("Seidensticker2017-", "", OBJ)) %>% 
  dplyr::mutate(label = gsub("Seidensticker2021-", "", label)) %>% 
  dplyr::mutate(path = paste0("output/gmm_masks/", OBJ, ".jpg")) %>%
  ggplot(aes(x = PC1, y = PC2, label = OBJ, image = path)) + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_image(size=.1) + 
  #geom_point() + 
  coord_equal() + 
  theme_bw()
ggsave("output/Fig_GMM_PCA.jpg", width = 8, height = 6)
#ggsave("Fig_GMM_PCA.pdf", width = 8, height = 6)
