# GMM prep

source("scripts/head.R")

src <- data.table::fread("https://raw.githubusercontent.com/dirkseidensticker/aSCAC/master/potterydrawings.csv")

# retrieve orig measurments from nwCongoDB
library(RSQLite)
temp <- tempfile()
utils::download.file(
  "https://github.com/dirkseidensticker/nwCongo/raw/master/data/base/nwCongoDB.sqlite", 
  temp, 
  mode = "wb", 
  quiet = TRUE)
con <- DBI::dbConnect(SQLite(), dbname = temp)

# get data from nwCongoDB
meas <- DBI::dbGetQuery(
  con,
  "SELECT 
    objID, muendungsD, maxD, minD 
  FROM t_Obj") %>%
  dplyr::mutate_if(is.character, as.numeric) %>%
  dplyr::mutate_all(~dplyr::na_if(., 0))


#BiocManager::install("EBImage")
library(EBImage)
#devtools::install_github('ISAAKiel/shapAAR')
library(shapAAR)

# Create Environment for multi-core computation
library(doParallel)
library(parallel)
cl <- makeCluster(detectCores() - 1)
print(cl)
registerDoParallel(cl)
on.exit(stopCluster(cl))

# list of new files:
f <- setdiff(
  tools::file_path_sans_ext(list.files("input/gmm_input/")),
  tools::file_path_sans_ext(list.files("output/gmm_masks/"))) # no need to redo already done files 

# combine measurements from DB with list of files
meas <- src %>% 
  dplyr::filter(SOURCE %in% f) %>% 
  dplyr::rename(objID = dbID) %>%
  dplyr::select(-c(grep("muendungsD", colnames(.)):grep("bodenD", colnames(.)))) %>%
  dplyr::left_join(meas, by = "objID")


foreach(
  i = 1:length(f),
  .packages = c(
    'shapAAR',
    'EBImage')) %dopar% {
      
      #meas.i <- meas[SOURCE == f[i],]
      
      # see https://github.com/ISAAKiel/shapAAR/tree/master/vignettes
      img <- EBImage::readImage(paste0("input/gmm_input/", f[i], ".png"))
      
      img <- getFrame(img, 1) # remove possible additional frames/layers
      img <- channel(img, "grey")
      
      # enhancing ----
      img <- EBImage::normalize(img)^5
      
      img <- EBImage::resize(img, dim(img)[1]/10) # dynamic resize of 10% size, keeps relative proportions?
      
      #img <- EBImage::resize(img, w = max(c(meas.i$muendungsD, meas.i$maxD, meas.i$minD), na.rm = T) * 10)
      
      img <- shapAAR::add_canvas(img,10,10,center = T)
      img <- img > EBImage::otsu(img)
      img <- EBImage::gblur(img,1)
      img <- round(1 + img * 255)
      img <- img - mean(img)
      g <- shapAAR::stopping_fun(img)
      
      # feature extraction ----
      v <- 1 # velocity of the active contour
      dt <- 1 # time interval that is used for the development of the contour
      phi <- default_phi(g) # starting contour = a 5 px border around the image
      n_iter <- 1000 # maximum iteration length of the algorithm
      buffer <- 2 # buffer to avoid that the contour would penetrate dashed lines
      
      phi_out <- active_contour(phi, 
                                g, 
                                n_iter = n_iter, 
                                v = v, 
                                dt = dt, 
                                show = FALSE, 
                                buffer = buffer)
      
      
      inner <- EBImage::Image(phi_out<=0)
      inner <- EBImage::fillHull(inner)
      for (j in 1:buffer) {
        inner <- EBImage::erode(inner)
      } # shrink the resulting contour by 3 pixel
      labelled_img <- EBImage::bwlabel(inner)
      features <- EBImage::computeFeatures.shape(labelled_img)
      all_objects <- 1:max(labelled_img)
      biggest_object <- which.max(features[,1])
      img_biggest_only <- EBImage::rmObjects(labelled_img,
                                             all_objects[all_objects != biggest_object])
      
      # feature optimisation ----
      #img_biggest_only <- shapAAR::img_crop_backgroud(img_biggest_only)
      
      # rectifying etc. ----
      fg_points <- which(img_biggest_only!=0,arr.ind = T)
      minbbox <- shapAAR::getMinBBox(fg_points)
      moment <- 90 - minbbox$angle
      if (abs(moment)>45) {moment = moment - sign(moment)*90}
      img_rect <- EBImage::rotate(img_biggest_only, moment, bg.col="black")
      img_rect <- img_rect > EBImage::otsu(img_rect)
      #img_rect <- shapAAR::img_crop_backgroud(img_rect)
      
      # save as jpg for Momocs ----
      writeImage(max(img_rect) - img_rect,
                 paste0("output/gmm_masks/", f[i], ".jpg"))
                 #paste0("output/gmm_masks/", gsub("\\..*","",f[i]), ".jpg"))
    }
