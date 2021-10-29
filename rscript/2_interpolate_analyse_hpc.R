library(Rmpi)
library(snow)

setwd('/home/tcrnbgh/Scratch/quarry_data/quarry_hpc')
sink('./2_interpolate_analyse_hpc_sinkout.txt')

# # check if packages need installing
# print('checking packages...')
# .libPaths(c('/home/tcrnbgh/R/x86_64-pc-linux-gnu-library/4.1',
#             .libPaths()))
# list.of.packages <- c("tuneRanger","gstat","reshape2",
#                       "e1071","caret","randomForest","stringr",
#                       "stars","sf","dplyr","gdalUtils",
#                       "raster","automap","fields","interp",
#                       "mgcv","purrr","furrr","doParallel",
#                       "future.apply", "snow")
# new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# if(length(new.packages)) {
#   print('following packages need installing...')
#   print(new.packages)
#   print('installing...')
#   install.packages(new.packages,
#                    lib="/home/tcrnbgh/R/x86_64-pc-linux-gnu-library/4.1")
#   print('done!')
# } else print('no new packages need installing!')
# print('done!')

# source functions
# print('loading functions...')
# 
# source('rscript/general_functions.R')
# print('done!')

# prepare data --------------------------------
print('loading prepped data...')
f <- 'data/prepData_alllocs_maxdiff01_smpper0.RDS'
prepData <- readRDS(f)
print('done!')
print('truncating prepData...')
prepDataTrunc <- prepData[1:70]
print('done!')

# interpolation run --------------------------------------------
print('running interpolations...')
st <- Sys.time()

#!! if offSet = T the use pd$tiles$pol !!
print('getting MPI cluster...')
cl <- snow::getMPIcluster()
print('done!')
# Display info about each process in the cluster
print(clusterCall(cl, function() Sys.info()))

print('running interpolations...')

datOut <- snow::clusterApply(cl, prepDataTrunc, function(pd) {
  
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
  setwd('/home/tcrnbgh/Scratch/quarry_data/quarry_hpc')
  print('temp dir...')
  print(tempdir())
  
  userDataDir <- '/home/tcrnbgh/Scratch/quarry_data'
  grassGISDBASE <- paste0(userDataDir,'/grassdb')
  f <- 'data/prepData_alllocs_maxdiff01_smpper0.RDS'
  sessionTag <- str_replace(basename(f),'.RDS','')
  sessionTag <- paste0(sessionTag,'')
  
  sink(paste0('logs/2_interpolate_analyse_hpc_sinkout_site',pd$pol$fid,'.txt'))
  
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
  
    # interpolation function -----
    
    maskPoly = pd$pol
    paramData = cvGrids
    gLoc = grassGISDBASE
    outputDir = '/home/tcrnbgh/Scratch/quarry_data/data_output'
    testCV = T # = T for test run
    outputTag = sessionTag
    intMethods=c(
     'rfsp',
     'nn','idw','ok','tin',
     'gfilter',
     'gspline'
    )

    trainingData <- pd$foldA$train
    testData <- pd$foldA$all
    # object for storing run times
    intTimes <- list()
    # object for storing interpolated rasters
    rasterlist <- list()
    # interpolation models ----
    
    # no parameters needed to be trialled for RF and TIN
    # random forest model, uses buffer distances----
    fit_RF <- function(trainingSp, testSp, tag='',subSample=NULL,
                       classInt=0.01) {
      # save crs from training data
      spCrs <- crs(trainingSp)
      
      # sample hole points
      if (!is.null(subSample)) {
        trainingSp <- trainingSp[sample(nrow(trainingSp),subSample),]
      }
      
      intTemplateCopy <- testSp
      # convert full a_ras to grid of locations to predict elev distances against
      gridPix <- as(testSp, "SpatialPixelsDataFrame")
      
      # split elevation into classes to reduce comp time
      classesElev <- cut(trainingSp$elev,
                         breaks=seq(min(trainingSp$elev),
                                    max(trainingSp$elev),
                                    by=classInt))
      print(paste0('splitting elevation values by ',
                   classInt, ' into ',length(levels(classesElev)),
                   ' classes...'))
      print('calculating buffer distances...')
      maxSize <- max(ncol(testData$ras),nrow(testData$ras))^2
      
      # grid.dist0 <- GSIF::buffer.dist(trainingSp["elev"], 
      #                                 gridPix[1], 
      #                                 # classes=as.factor(1:length(unique(trainingSp$elev))), # all elevs
      #                                 classes=classesElev,# groups elev values
      #                                 width=maxSize # max search radius
      # )
      # 
      grid.dist0 <- buffer.dist2(trainingSp["elev"],
                                 gridPix[1],
                                 # classes=as.factor(1:length(unique(trainingSp$elev))), # all elevs
                                 classes=classesElev,# groups elev values
                                 width=maxSize # max search radius
      )
      
      buffDists <- cbind(grid.dist0@data,gridPix@data)
      grid.dist0@data <- buffDists
      
      ov.elev <- over(trainingSp["elev"], grid.dist0)
      
      
      # find best mtry value
      print('tuning Random Forest to find best "mtry" value...')
      x <- ov.elev[,setdiff(names(ov.elev),'elev')]
      y <- ov.elev[,'elev']
      dn0 <- paste(setdiff(names(ov.elev),c('elev','slope','aspect')),
                   collapse="+")
      fm0 <- as.formula(paste("elev ~ ", dn0))
      
      # tune and construct RF
      m.elev <- tuneRanger::tuneMtryFast(fm0, 
                                         data = ov.elev,
                                         num.treesTry = 250,
                                         stepFactor = 1.5, 
                                         improve = 1e-5,
                                         trace = TRUE, 
                                         plot = FALSE,
                                         doBest = T,
                                         quantreg=TRUE)
      
      elev.rfd <- predict(m.elev, grid.dist0@data, type="quantiles")$predictions
      
      intTemplateCopy$pred_elev <- elev.rfd[,2]
      
      outList <- list(rfSp = rasterFromXYZ(cbind(intTemplateCopy@coords,
                                                 intTemplateCopy$pred_elev),
                                           crs=spCrs))
      
      return(outList)
    }
    
    if ('rfsp' %in% intMethods) {
      st <- Sys.time()
      interp_RF <- fit_RF(trainingSp = trainingData$sp,
                          testSp = testData$sp,
                          tag='classed')
      
      intTimes$RF <- Sys.time()-st
      rasterlist$`Random Forest SP` <- list(interp_RF$rfSp)
      cat('completed rfSp', file=paste0('rfSp_',pd$pol$fid,'.txt'))
    }
    
    if ('tin' %in% intMethods) {
      # Triangular Irregular Surface
      st <- Sys.time()
      fit_TIN <- interp::interp( # using {interp}
        x = trainingData$df$X,     # the function actually accepts coordinate vectors
        y = trainingData$df$Y,
        z = trainingData$df$elev,
        xo = testData$df$X, # here we already define the target grid
        yo = testData$df$Y,
        output = "points") %>%
        bind_cols()
      
      interp_TIN <- raster::rasterFromXYZ(fit_TIN,
                                          crs = crs(trainingData$sf))
      intTimes$TIN <- Sys.time()-st
      
      rasterlist$`Triangular Irregular Surface` <- list(interp_TIN)
      cat('completed TIN', file=paste0('TIN_',pd$pol$fid,'.txt'))
    }
    
    # parameterised interpolation models
    # process CV grids to remove bad combinations (nmin vals > nmax vals)
    # and add run number, add polygon id,
    # test routine ensures parameter with greatest impact on computation time is
    # evenly sampled.
    # df <- paramData$nn
    # paramData <- cvGrids
    
    paramData.c <- 
      paramData %>% 
      map2(.y = names(paramData), 
           .f = function(df, y) {
             if (any(str_detect(names(df),'nmaxVals')) && 
                 any(str_detect(names(df),'nminVals'))) {
               df <- df %>% filter(nminVals < nmaxVals) }
             df <- df %>% 
               mutate(run_no = 1:nrow(.),
                      intpol_fid = maskPoly$fid,
                      int_method = y)
             if (testCV) {
               if (y == 'nn' | y == 'idw' | y == 'ok') {
                 df <- rbind(df[which.max(df$nmaxVals),],
                             df[which.min(df$nmaxVals),],
                             df[which.median(df$nmaxVals),]
                 ) }
               if (y == 'gspline') {
                 df <- rbind(df[which.max(df$nminVals),],
                             df[which.min(df$nminVals),],
                             df[which.median(df$nminVals),]
                 ) }
               if (y == 'gbicubic') {
                 df <- rbind(df[which.max(df$stepVals),],
                             df[which.min(df$stepVals),],
                             df[which.median(df$stepVals),]
                 ) }
               if (y == 'gfilter') {
                 df <- rbind(df[which.max(df$radVals),],
                             df[which.min(df$radVal),],
                             df[which.median(df$radVals),]
                 ) }
             }
             return(df)
           })
    
    if ('nn' %in% intMethods) {
      interp_NNs <- list()
      for (y in 1:length(paramData.c$nn$run_no)) {
        x <- paramData.c$nn$run_no[y]
        pdata <- paramData.c$nn %>% 
          filter(run_no == x)
        print(pdata)
        # Nearest neighbour
        st <- Sys.time()
        fit_NN <- gstat::gstat( # using package {gstat} 
          formula = elev ~ 1,    
          data = trainingData$sp,
          nmax = pdata$nmaxVals,
          nmin = pdata$nminVals,# Number of neighboring observations used for the fit
        )
        interp_NNs[[y]] <- raster::interpolate(testData$ras[[1]], fit_NN)
        tdiff <- Sys.time()-st
        t <- list(val = tdiff,
                  unit_chr = units(tdiff))
        intTimes$NN[[y]] <- t
      }
      rasterlist$`Nearest Neighbor` <- interp_NNs
      cat('completed NN', file=paste0('NN_',pd$pol$fid,'.txt'))
    }
    
    if ('idw' %in% intMethods) { 
      interp_IDWs <- list()
      for (y in 1:length(paramData.c$idw$run_no)) {
        x <- paramData.c$idw$run_no[y]
        pdata <- paramData.c$idw %>% 
          filter(run_no == x)
        
        # Inverse Distance Weighting
        st <- Sys.time()
        fit_IDW <- gstat::gstat( # The setup here is quite similar to NN
          formula = elev ~ 1,
          data = trainingData$sp,
          nmax = pdata$nmaxVals,
          nmin = pdata$nminVals,
          set = list(idp = pdata$idpVals)) # inverse distance power
        
        interp_IDWs[[y]] <- raster::interpolate(testData$ras[[1]], fit_IDW)
        tdiff <- Sys.time()-st
        t <- list(val = tdiff,
                  unit_chr = units(tdiff))
        intTimes$IDW[[y]] <- t
      }
      rasterlist$`Inverse Distance Weighted` <- interp_IDWs
      cat('completed IDW', file=paste0('IDW_',pd$pol$fid,'.txt'))
    }
    
    if ('ok' %in% intMethods) {
      ## ordinary kriging
      v_mod_ok <- automap::autofitVariogram(elev~1,input_data =
                                              trainingData$sp)
      
      # not happy running with future_lapply!
      interp_OKs <- list()
      for (y in 1:length(paramData.c$ok$run_no)) {
        x <- paramData.c$ok$run_no[y]
        pdata <- paramData.c$ok %>% 
          filter(run_no == x)
        # print(OK_grid[x,])
        st <- Sys.time()
        fit_OK = gstat(formula = elev ~ 1, model = v_mod_ok$var_model, 
                       data = trainingData$sp,
                       nmax=pdata$nmaxVals,
                       nmin=pdata$nminVals)
        interp_OKs[[y]] <- raster::interpolate(testData$ras[[1]], fit_OK)
        tdiff <- Sys.time()-st
        t <- list(val = tdiff,
                  unit_chr = units(tdiff))
        intTimes$OK[[y]] <- t
      }
      rasterlist$`Ordinary Kriging` <- interp_OKs
      cat('completed OK', file=paste0('OK_',pd$pol$fid,'.txt'))
    }
    
    # grass params
    # generate unique grass location name
    gnLoc <- paste0(gLoc,'/grass_pol',pd$pol$fid)
    if (dir.exists(gnLoc)) unlink(gnLoc,recursive=T)
    gdb <- paste0(gnLoc,'/PERMANENT')
    
    # grass-based splines model ----
    # training.sf <- trainingData$sf
    # test.r <- testData$ras
    
    if ('gspline' %in% intMethods) {
      print('attempting gspline vrts interpolation...')
      # set up gdb
      # make location with grass call
      system(paste0('grass -c /home/tcrnbgh/quarry_hpc/vector/init_vector.gpkg ',gnLoc,' -e'))
      system(paste0(
        'grass ',gdb,' --exec g.proj datum=osgb36 -c'
      ))
      print(gdb)
      # write points for gspline
      if (file.exists(paste0('vector/intout_',maskPoly$fid,'_training.gpkg'))) {
        file.remove(paste0('vector/intout_',maskPoly$fid,'_training.gpkg'))
      }
      st_crs(trainingData$sf) <- st_crs(27700)
      
      st_write(trainingData$sf[,'elev'], 
               paste0('vector/intout_',maskPoly$fid,'_training.gpkg'))
      vecLoc <- paste0(getwd(),'/vector/intout_',maskPoly$fid,'_training.gpkg')
      
      system(paste0(
        'grass ',gdb,' --exec v.in.ogr input=',vecLoc,' output=points -o --o'
      ))
      
      m <- mask(testData$ras[[1]],trainingData$ras[[1]],inverse=T)
      mLoc <- paste0(getwd(),'/raster/intout_',maskPoly$fid,'_ras_testmask.tif')
      
      m.new <- st_as_stars(m) %>% 
        st_set_crs(st_crs(27700)) %>% as(.,'Raster')
      print('writing mask data as tif...')
      writeRaster(m.new, mLoc,
                  overwrite=T)
      print('done!')
      
      
      system(paste0(
        'grass ',gdb,' --exec r.in.gdal input=',mLoc,' output=testMask -o --o'
      ))
      system(paste0(
        'grass ',gdb,' --exec g.region raster=testMask'
      ))
      # y <- 1
      interp_GSPLINEs <- list()
      for (y in 1:length(paramData.c$gspline$run_no)) {
        x <- paramData.c$gspline$run_no[y]
        pdata <- paramData.c$gspline %>% 
          filter(run_no == x)
        st <- Sys.time()
        system2('grass',
                paste(shQuote(gdb),
                      shQuote('--exec'),
                      shQuote(paste0(getwd(),'/python/GRASS_vrst.py')),
                      shQuote(pdata$smoothVals),
                      shQuote(pdata$tensionVals),
                      shQuote(pdata$nminVals),
                      shQuote(x),
                      shQuote(pdata$intpol_fid),
                      shQuote(paste0(getwd(),'/raster'))
                ),
                stderr = paste0(getwd(),'/logs/grass_vrst_errout_',
                                pdata$intpol_fid,'.txt')
        )
        print('finished executing grass call...')
        print('reading grass output raster...')
        if (file.exists(paste0('raster/gspline_int_intfid_',pdata$intpol_fid,
                               '_runnum_',x,'.tif'))) {
          
        r <- readAll(raster(paste0('raster/gspline_int_intfid_',pdata$intpol_fid,
                           '_runnum_',x,'.tif')))
        print('done!')
        interp_GSPLINEs[[y]] <- raster::merge(r,trainingData$ras[[1]])
        # file.remove(paste0('raster/gspline_int_intfid_',pdata$intpol_fid,
        #                    '_runnum_',x,'.tif'))
        tdiff <- Sys.time()-st  
        t <- list(val = tdiff,
                  unit_chr = units(tdiff))
        intTimes$GSPLINE[[y]] <- t
        } else { 
          print(paste0('failed interpolation...pol_fid: ',pd$pol$fid,' run no: ',
                       y))
          interp_GSPLINEs[[y]] <- NULL 
        }
      }
      rasterlist$`GRASS Regularized Splines Tension` <- interp_GSPLINEs
      cat('completed GSPLINE', file=paste0('GSPLINE_',pd$pol$fid,'.txt'))
      if (dir.exists(gnLoc)) unlink(gnLoc,recursive=T)
    }
    
    
    if ('gfilter' %in% intMethods) {
      print('attempting gfilter interpolation...')
      # make location with grass call
      system(paste0('grass -c /home/tcrnbgh/quarry_hpc/vector/init_vector.gpkg ',gnLoc,' -e'))
      system(paste0(
        'grass ',gdb,' --exec g.proj datum=osgb36 -c'
      ))
      print(gdb)
      
      
      # prepare rasters
      mLoc <- paste0(getwd(),'/raster/intout_',maskPoly$fid,'_ras_training.tif')
      # not sure why rasters were merged??
      # allRas <- merge(trainingData$ras[[1]],testData$ras[[1]]) # weird?
      # overwritten below:
      
      allRas <- trainingData$ras[[1]] %>% st_as_stars %>% 
        st_set_crs(st_crs(27700))
      
      print('writing training data as tif (stars method)...')
      write_stars(allRas, 
                  mLoc)
      print('done!')
      system(paste0(
        'grass ',gdb,' --exec r.in.gdal input=',mLoc,' output=allRas -o --o'
      ))
      system(paste0(
        'grass ',gdb,' --exec g.region raster=allRas'
      ))
      # x <- 1
      # y <- 1
      load(file='cvdev/gfilter.RDS')
      interp_GFILTERs <- list()
      for (y in 1:length(paramData.c$gfilter$run_no)) {
        x <- paramData.c$gfilter$run_no[y]
        
        pdata <- paramData.c$gfilter %>% 
          filter(run_no == x)
        
        # get filter text info
        pdata$filtText <- 
          gfilter.params[which(gfilter.params$filtcomb_number == pdata$filtVals),'col3']
        pdata$filtLength <-
          unlist(lapply(str_split(pdata$filtText,','),length))
        
        pdata$nradVals <- 
          ifelse(pdata$filtLength == 2, 
                 paste0(rep(pdata$radVals,2),collapse=','),
                 as.character(pdata$radVals))
        
        st <- Sys.time()
        system2('grass',
                paste(shQuote(gdb),
                      shQuote('--exec'),
                      shQuote(paste0(getwd(),'/python/GRASS_resampfilter.py')),
                      shQuote(pdata$nradVals),
                      shQuote(pdata$filtText),
                      shQuote(x),
                      shQuote(pdata$intpol_fid),
                      shQuote(paste0(getwd(),'/raster'))
                ),
                stderr = paste0(getwd(),'/logs/grass_resampfilter_errout_',
                                pdata$intpol_fid,'.txt')
        )
        print('finished executing grass call...')
        print('reading raster...')
        # for some reason raster() doesnt load data into R memory!
        if (file.exists(paste0('raster/gfilter_int_intfid_',pdata$intpol_fid,
                               '_runnum_',x,'.tif'))) {
          r <- readAll(raster(paste0('raster/gfilter_int_intfid_',pdata$intpol_fid,
                                     '_runnum_',x,'.tif')))
          print('done!')
          crs(r) <- crs(testData$ras[[1]])
          
          interp_GFILTERs[[y]] <- r
          # file.remove(paste0('raster/gfilter_int_intfid_',pdata$intpol_fid,
          #                    '_runnum_',x,'.tif'))
          tdiff <- Sys.time()-st
          t <- list(val = tdiff,
                    unit_chr = units(tdiff))
          intTimes$GFILTER[[y]] <- t
        } else { 
          print(paste0('failed interpolation...pol_fid: ',pd$pol$fid,' run no: ',
                       y))
          interp_GFILTERs[[y]] <- NULL 
          }
      }
      rasterlist$`GRASS Resampled Filter` <- interp_GFILTERs
      cat('completed GFILTER', file=paste0('GFILTER_',pd$pol$fid,'.txt'))
      if (dir.exists(gnLoc)) unlink(gnLoc,recursive=T)
    }
    
    # gen list
    rasterlist <- rasterlist %>% 
      map( ~map(.x, .f = function(r) {
        mask(crop(extend(r,testData$ras[[1]]),testData$ras[[1]]),
             testData$ras[[1]])
      }))
    cat('completed raster list', file=paste0('rasterlist_',pd$pol$fid,'.txt'))
    
    # output parameters as melted df
    paramsCv <- paramData.c %>% 
      map_df(~reshape2::melt(.x,
                             id.vars=c('run_no','intpol_fid',
                                       'int_method')))
    
    print(paste0('completed pol id...',maskPoly$fid))
    intA <- list(ras = rasterlist,
                 params.m = paramsCv,
                 intTimes = intTimes)
    
    dat <- compareInt(intRasters=intA,
                      foldedRas=pd$foldA,
                      tiledRas=pd$tiles)
    cat('completed compareInt', file=paste0('compareInt_',pd$pol$fid,'.txt'))
    
    fout <- paste0(outputDir,'/intdat_',outputTag,'_polfid',pd$pol$fid,'.RDS')
    save(dat,
         file=fout)
    print(paste0('saved file to... ',fout))
    intTimes
})
print('done!')
# Clean up the cluster and release the relevant resources.
stopCluster(cl)
sink()
mpi.quit()


# clusterEvalQ (cl, traceback())

# 
# system('cp /home/tcrnbgh/quarry_hpc/rscript/* /home/tcrnbgh/Scratch/quarry_data/quarry_hpc/rscript/')
# system('cp /home/tcrnbgh/quarry_hpc/python/* /home/tcrnbgh/Scratch/quarry_data/quarry_hpc/python/')
# source(paste0(getwd(),'/rscript/general_functions.R'))
# pd <- prepData[[1]]
# st <- Sys.time()
# intTimes <- interpolateRas(pd,
#                maskPoly = pd$pol,
#                paramData=cvGrids,
#                outputDir = '/home/tcrnbgh/Scratch/quarry_data/data_output',
#                gLoc = grassLocation,
#                testCV = F, # = T for test run
#                outputTag = sessionTag,
#                intMethods=c(
#                  'rfsp',
#                  'nn','idw','ok','tin',
#                  'gfilter',
#                  'gspline'
#                )
# )
# intTimeDiff <- Sys.time() - st
# 
# 
# # clear grass mapsets
# allMaps <- list.dirs(substr(grassLocation,1,nchar(grassLocation)-1),
#                      recursive = F, full.names = T)
# dirRem <- setdiff(allMaps,paste0(grassLocation,"PERMANENT"))
# map(dirRem,~unlink(.x,recursive=T))
# # clear tmp rasters
# f <- list.files('raster',full.names = T)
# file.remove(f)


# # # bicubic cant be run in parallel so separate function is used here
# biOut <- lapply(prepData,interpolateRasBicubic,
#                 cvg=cvGrids,
#                 testCV=F, # check testCV = F for full run
#                 tag=sessionTag)
# 

# 
# bindCubic(sessionTag)
# print(Sys.time() - st)
