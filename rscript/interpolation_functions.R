library(tuneRanger)
library(gstat)
library(reshape2)
library(e1071)
library(caret)
library(randomForest)
library(stringr)
library(stars)
library(sf)
library(sp)
library(dplyr)
library(gdalUtils)
library(raster)
library(automap)
library(fields)
library(interp)
library(mgcv) 
library(purrr)
# new CV grids ----------------------------------------------------
# prep parameters for grass resamp.filter 
print('generating cv params...')
load(file='cvdev/gfilter.RDS')

# revised cv parameters based on results of exploratory values above
cvGrids <- 
  list(
    nn = data.frame(expand.grid(
      nmaxVals = seq(1,100,by=2),
      nminVals = 1
    )),
    idw = data.frame(expand.grid(
      nmaxVals = seq(2,100,by=2),
      nminVals = 1,
      idpVals = seq(0.2,10,by=0.2)
    )),
    ok = data.frame(expand.grid(
      nmaxVals = seq(1,100,by=2),
      nminVals = 1
    )),
    gspline = data.frame(expand.grid(
      tensionVals = c(seq(0.0001,0.1,by=0.002),
                      seq(0.15,0.7,by=0.05)),
      smoothVals = seq(2,30,by=5),
      nminVals = seq(30,150,by=20)
    )),
    gbicubic = data.frame(expand.grid(
      stepVals = seq(1,20,by=1),
      lamVals = c(seq(0.0001,0.1,by=0.005),seq(0.1,1.4,by=0.1))
    )),
    gfilter = data.frame(expand.grid(
      radVals = seq(5,50,by=5),
      filtVals = gfilter.params$filtcomb_number
    ))
  )
print('done!')

# end new CV grids ------------------------------------------------

which.median <- function(x) which.min(abs(x - median(x)))

# useful print paste function
pp <- function(...) {
  print(paste0(list(...),collapse=''))
}

'%!in%' <- function(x,y)!('%in%'(x,y))

# intRasters <- intA
# foldedRas <- pd$foldA
# tiledRas <- pd$tiles

compareInt <- function(intRasters, # list of interpolated rasters
                       foldedRas, # the test/training rasters (raster A, latest)
                       tiledRas # another raster to compare int surfaces with (raster B, earliest)
) { 
  
  # the earlier survey to compare against
  compareRas <- tiledRas$b
  # a polygon to exclude from testErr calculations
  maskPoly <- tiledRas$pol
  
  frtest <- as(foldedRas$test,'Raster')
  
  foldedRas$test <- frtest[[1]]
  
  if (class(compareRas)=='stars') {
    compareRas.r <- as(compareRas, 'Raster')
    compareRas.r <- compareRas.r[[1]]
  }
  
  # this ensures the extent of the comparison raster (older survey)
  # is the same as the interpolated surfaces. Errors can occur when
  # earlier survey did not quite cover the same cells as the later one.
  compareRas.r <- mask(crop(extend(compareRas.r,foldedRas$all$ras[[1]]),
                            foldedRas$all$ras[[1]]),
                       foldedRas$all$ras[[1]])
  
  # training error is difference between original training data and modeled training data
  # test error is difference between original test data and modeled test data
  
  # test data below is essentially the original raster 
  # ep <- pd$pol
  # fr <- foldedRas
  # cr <- compareRas.r
  # raslist <- intRasters$ras[[1]]
  compareEach <- function(raslist,fr,cr,ep=NULL) {
    
    # raslist <- intRasters$ras$`GRASS Regularized Splines Tension`
    # interpolated <- raslist[[1]]
    compFunction <- function(interpolated) {
      
      trainingErr.r <- mask(interpolated,fr$train$ras[[1]]) - 
        fr$train$ras[[1]]
      
      testErr.r <- mask(interpolated,fr$test) - fr$test
      compareDiff <- interpolated - cr
      
      out <- list(
        trainingErr.r = trainingErr.r,
        testErr.r = testErr.r,
        compareDiff = compareDiff)
      
      if (!is.null(ep)) {
        ep.r <- fasterize::fasterize(ep, fr$test)
        out$testErr.ex.r <- mask(interpolated,fr$test) - 
          mask(fr$test,ep.r,inverse=T)
        out$testErr.inc.r <- mask(interpolated,fr$test) - 
          mask(fr$test,ep.r)
        out$compareDiff.inc.r <- interpolated - mask(cr,ep.r)
        out$compareDiff.ex.r <- interpolated - mask(cr,ep.r,inverse=T)
      }
      
      return(out)
    }
    
    updatedRas <- raslist %>% map(~compFunction(.x))
    # raslist$ras.comp <- updatedRas
    # return(raslist)
  }
  
  diffMaps <- intRasters$ras %>% 
    map(~compareEach(.x, 
                     fr = foldedRas,
                     cr = compareRas.r,
                     ep = maskPoly))
  
  # diffMaps$`Nearest Neighbor`$compare$compareDiff
  # plot(stack(diffMaps$`Triangular Irregular Surface`))
  # r <- diffMaps$`Ordinary Kriging`$trainingErr.r
  calcRMSEfromRas <- function(r) sqrt(cellStats(r^2,mean))
  
  diffRMSEs <- lapply(names(diffMaps), function(x) {
    diffRMSEs <- diffMaps[[x]] %>% map_df(map, calcRMSEfromRas) %>% 
      mutate(int_method = x,
             run_no = 1:nrow(.))
  }) %>% bind_rows()
  
  # add pol fid if pol provided
  if (!is.null(maskPoly)) diffRMSEs <- diffRMSEs %>% 
    mutate(intpol_fid = maskPoly$fid)
  
  return(list(orig.maps = intRasters,
              # diff.maps = diffMaps,
              tiles = tiledRas,
              rmses = diffRMSEs))
}


buffer.dist2 <- function(observations, predictionDomain, classes, width, ...) {
  if(missing(width)){ width <- sqrt(areaSpatialGrid(predictionDomain)) }
  if(!length(classes)==length(observations)){ stop("Length of 'observations' and 'classes' does not match.") }
  ## remove classes without any points:
  xg = summary(classes, maxsum=length(levels(classes)))
  selg.levs = attr(xg, "names")[xg > 0]
  
  if(length(selg.levs)<length(levels(classes))){
    fclasses <- as.factor(classes)
    fclasses[which(!fclasses %in% selg.levs)] <- NA
    classes <- droplevels(fclasses)
  }
  
  ## derive buffer distances
  s <- list(NULL)
  
  for(i in 1:length(levels(classes))) {
    if (nrow(observations[which(classes==levels(classes)[i]),1]) > 0) {
      
      s[[i]] <- raster::distance(
        rasterize(
          observations[which(classes==levels(classes)[i]),1]@coords, y=raster(predictionDomain)),
        width=width)
    } else print(i)
  }
  s <- s[sapply(s, function(x){!is.null(x)})]
  s <- brick(s)
  s <- as(s, "SpatialPixelsDataFrame")
  s <- s[predictionDomain@grid.index,]
  return(s)
}