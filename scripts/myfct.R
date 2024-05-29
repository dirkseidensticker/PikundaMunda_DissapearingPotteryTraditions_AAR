rcarbonsum <- function(a, 
                       oxcalnorm = FALSE){
  
  cal = rcarbon::calibrate(x = a$C14AGE,
                           errors = a$C14STD,
                           calCurves = 'intcal20', 
                           ncores = ncores, 
                           normalised = FALSE) #running calibration over 3 cores
  
  spd <- rcarbon::spd(cal,
                      timeRange=c(4000,0), 
                      spdnormalised = TRUE)
  
  spd <- as.data.frame(spd[2])
  
  if(oxcalnorm) {
    # raise to max() == 1 like OxCal does!
    spd$grid.PrDens <- spd$grid.PrDens/max(spd$grid.PrDens, na.rm = TRUE)
  }
  
  spd <- spd[spd$grid.PrDens != 0,] # remove 0 values
  spd$grid.calBP <- 1950 - spd$grid.calBP
  
  med <- list()
  for(k in 1:length(cal$grids)){
    m <- 1950 - median(cal$grids[[k]]$calBP, na.rm = T)
    med[k] <- m
  }
  median <- do.call(rbind, med)
  
  start <- min(spd$grid.calBP)
  
  output <- list(dat = spd, 
                 median = median, 
                 start = start)
}
