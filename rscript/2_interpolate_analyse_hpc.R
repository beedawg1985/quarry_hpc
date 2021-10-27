library(Rmpi)
library(snow)

setwd('/home/tcrnbgh/Scratch/quarry_data/quarry_hpc')
sink('./2_interpolate_analyse_hpc_sinkout.txt')

# check if packages need installing
print('checking packages...')
.libPaths(c('/home/tcrnbgh/R/x86_64-pc-linux-gnu-library/4.1',
            .libPaths()))
list.of.packages <- c("tuneRanger","gstat","reshape2",
                      "e1071","caret","randomForest","stringr",
                      "stars","sf","dplyr","gdalUtils",
                      "raster","automap","fields","interp",
                      "mgcv","purrr","furrr","doParallel",
                      "future.apply", "snow")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages,
                                          lib="/home/tcrnbgh/R/x86_64-pc-linux-gnu-library/4.1")

print('done!')
# set envs
userDataDir <<- '/home/tcrnbgh/Scratch/quarry_data'
grassMapset <<- paste0(userDataDir,'/grassdb/quarry/PERMANENT/')
grassLocation <<- paste0(userDataDir,'/grassdb/quarry/')

# source functions
print('loading functions...')

source('rscript/general_functions.R')
print('done!')

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

# prepared data --------------------------------
print('loading prepped data...')
f <- 'data/prepData_alllocs_maxdiff01_smpper0.RDS'
prepData <- readRDS(f)
print('done!')
# interpolation run --------------------------------------------

# these calls establish workers for parallel processing and retrieve the id of
# each worker. The id of the workers is tallied with a pre-existing GRASS mapset
# id, so each worker simultanesouly processes data in a single mapset to avoid 
# conflicts.
# run
sessionTag <- str_replace(basename(f),'.RDS','')
sessionTag <- paste0(sessionTag,'')
print('running interpolations...')
st <- Sys.time()

#!! if offSet = T the use pd$tiles$pol !!
print('getting MPI cluster...')

cl <- snow::getMPIcluster()

# Display info about each process in the cluster
print(clusterCall(cl, function() Sys.info()))

# if (is.null(cl)) { 
#   print('getMPIcluster() failed...')
#   print('trying to make cluster...')
#   cl <- 
#     snow::makeMPIcluster(count=35,
#                          outfile=paste0(getwd(),'/logs/cluster_out.txt'))
#   cl <- makeCluster(mc <- getOption("cl.cores", 35),
#                     outfile=paste0(getwd(),'/logs/cluster_out.txt'))
#   }
print('done!')

print('sourcing functions and exporting vars on cluster...')
# snow::clusterEvalQ(cl, source('/home/tcrnbgh/Scratch/quarry_data/quarry_hpc/rscript/general_functions.R'))
snow::clusterExport(cl, list=c('cvGrids','grassLocation','sessionTag',
                               'userDataDir','grassMapset'))
print('done!')

print('running interpolations...')
prepDataTrunc <- prepData[1:70]

datOut <- snow::clusterApply(cl, prepDataTrunc, function(pd) {
  
  interpolateRas <- function(pd,
                             maskPoly, 
                             paramData,
                             testCV = T,
                             outputDir = '/media/mal/working_files/quarry/',
                             outputTag='testnew',
                             gLoc='/home/barneyharris/user_quarry_data/grassdb/quarry/',
                             # gMap='q1',
                             intMethods=c('nn','idw','ok',
                                          'rfsp','tin','gspline',
                                          'gbicubic','gfilter') ) {
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
    }
    
    # grass params
    gdb <- paste0(gLoc,Sys.getpid(),'/')
    
    # @@@@ needed to prevent freezing on run 279 ?!? !!!!
    # paramData.c$gspline <- paramData.c$gspline %>%
    #   slice(10:nrow(.))
    # grass-based splines model ----
    # training.sf <- trainingData$sf
    # test.r <- testData$ras
    if ('gspline' %in% intMethods) { 
      # write points for gspline
      st_write(trainingData$sf[,'elev'], 
               paste0('vector/intout_',maskPoly$fid,'_training.gpkg'),
               delete_dsn=T)
      vecLoc <- paste0(getwd(),'/vector/intout_',maskPoly$fid,'_training.gpkg')
      
      system(paste0(
        'grass ',gdb,' --exec g.proj datum=osgb36 -c'
      ))
      
      if (dir.exists(gdb)) unlink(gdb,recursive=T)
      
      system(paste0(
        'grass -c ',gdb,' --exec v.in.ogr input=',vecLoc,' output=points -o --o'
      ))
      
      m <- mask(testData$ras[[1]],trainingData$ras[[1]],inverse=T)
      mLoc <- paste0(getwd(),'/raster/intout_',maskPoly$fid,'_ras_testmask.tif')
      writeRaster(m, mLoc,
                  overwrite=T)
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
        print('reading grass output raster...')
        r <- raster(paste0('raster/gspline_int_intfid_',pdata$intpol_fid,
                           '_runnum_',x,'.tif'))
        interp_GSPLINEs[[y]] <- raster::merge(r,trainingData$ras[[1]])
        file.remove(paste0('raster/gspline_int_intfid_',pdata$intpol_fid,
                           '_runnum_',x,'.tif'))
        tdiff <- Sys.time()-st
        t <- list(val = tdiff,
                  unit_chr = units(tdiff))
        intTimes$GSPLINE[[y]] <- t
      }
      rasterlist$`GRASS Regularized Splines Tension` <- interp_GSPLINEs
    }
    
    
    if ('gfilter' %in% intMethods) {
      # prepare rasters
      mLoc <- paste0(getwd(),'/raster/intout_',maskPoly$fid,'_ras_training.tif')
      # not sure why rasters merged??
      allRas <- merge(trainingData$ras[[1]],testData$ras[[1]]) # weird?
      # overwritten below:
      allRas <- trainingData$ras[[1]]
      
      writeRaster(allRas, 
                  mLoc,
                  overwrite=T)
      
      if (dir.exists(gdb)) unlink(gdb,recursive=T)
      
      system(paste0(
        'grass -c ',gdb,' --exec r.in.gdal input=',mLoc,' output=allRas -o --o'
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
        r <- readAll(raster(paste0('raster/gfilter_int_intfid_',pdata$intpol_fid,
                                   '_runnum_',x,'.tif')))
        crs(r) <- crs(testData$ras[[1]])
        
        interp_GFILTERs[[y]] <- r
        file.remove(paste0('raster/gfilter_int_intfid_',pdata$intpol_fid,
                           '_runnum_',x,'.tif'))
        tdiff <- Sys.time()-st
        t <- list(val = tdiff,
                  unit_chr = units(tdiff))
        intTimes$GFILTER[[y]] <- t
      }
      rasterlist$`GRASS Resampled Filter` <- interp_GFILTERs
    }
    
    # gen list
    rasterlist <- rasterlist %>% 
      map( ~map(.x, .f = function(r) {
        mask(crop(extend(r,testData$ras[[1]]),testData$ras[[1]]),
             testData$ras[[1]])
      }))
    
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
    fout <- paste0(outputDir,'/intdat_',outputTag,'_polfid',pd$pol$fid,'.RDS')
    save(dat,
         file=fout)
    print(paste0('saved file to... ',fout,
                 outputDir,'/intdat_',outputTag,'_polfid',pd$pol$fid,'.RDS'))
    return(intTimes)
  }
  
  interpolateRas(pd,
                 maskPoly = pd$pol,
                 paramData=cvGrids,
                 gLoc = grassLocation,
                 outputDir = '/home/tcrnbgh/Scratch/quarry_data/data_output',
                 testCV = T, # = T for test run
                 outputTag = sessionTag,
                 intMethods=c(
                   'rfsp',
                   'nn','idw','ok','tin',
                   'gfilter',
                   'gspline'
                 )
                 )
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
