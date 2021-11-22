library(Rmpi)
library(snow)

library(raster)
library(stars)
library(sf)
library(dplyr)
library(purrr)
setwd('/home/tcrnbgh/Scratch/quarry_data/quarry_hpc')
sink('./2_interpolate_analyse_hpc_sinkout.txt')

Sys.setenv(TMPDIR='/home/tcrnbgh/Scratch/tmp')
print('main R tempdir is...')
print(tempdir())

# prepare data --------------------------------
print('loading prepped data...')
prepData <- readRDS('data/prepData_alllocs_norm_maxdiff01_smpper0.RDS')

jackKnife <- function(pd) {
  newPd <- 1:nrow(pd$foldA$train$sf) %>%                          # check this!
  # newPd <- 1:10 %>% 
    map(.f = function(x) {
      newTrain <- list()
      pd.sf <- pd$foldA$train$sf %>% st_sf
      newTrain$sf <- pd.sf[setdiff(1:nrow(pd$foldA$train$sf),x),]
      newTrain$df <- cbind(st_coordinates(newTrain$sf),
                           st_drop_geometry(newTrain$sf))
      newTrain$sp <- as_Spatial(newTrain$sf)
      crs(newTrain$sp) <- crs(pd$foldA$train$sp)
      newTrain$ras <- rasterFromXYZ(newTrain$df)
      crs(newTrain$ras) <- crs(newTrain$sp)
      pd$foldA$train <- newTrain
      
      newTestRas <- mask(pd$foldA$train$ras,newTrain$ras,
                         inverse=T)
      crs(newTestRas) <- crs(newTrain$sp)
      pd$foldA$test <- st_as_stars(newTestRas)
      
      pd$pol$fid <- paste0(pd$pol$fid,'.',x)
      crs(pd$foldA$all$sp) <- crs(newTrain$sp)
      pd
    })
}

prepData <- prepData[1:5] %>% map(jackKnife) %>% 
  flatten()

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

datOut <- snow::clusterApply(cl, prepData, function(pd) {
  
  setwd('/home/tcrnbgh/Scratch/quarry_data/quarry_hpc')
  source('rscript/interpolation_functions.R')
  print('temp dir...')
  print(tempdir())
  Sys.setenv(TMPDIR='/home/tcrnbgh/Scratch/tmp')
  
  userDataDir <- '/home/tcrnbgh/Scratch/quarry_data'
  grassGISDBASE <- paste0(userDataDir,'/grassdb')

  sessionTag <- 'prepData_alllocs_norm_maxdiff01_jackknife'        # check this !
  
  sink(paste0('logs/2_interpolate_analyse_hpc_sinkout_site',
              pd$pol$fid,'.txt'))
  
  cvGrids <- loadCV()
  
  interpolateRas(pd,
    paramData = cvGrids,
    gLoc = grassGISDBASE,
    outputDir = '/home/tcrnbgh/Scratch/quarry_data/data_output',
    testCV = F, # = T for test run                                # check this !
    outputTag = sessionTag,
    intMethods=c(
      # 'rfsp','tin',
      'nn','idw','ok',
      'gbicubic',
      'gspline'
      )
    )
  
})
save(datOut,
     file=paste0('/home/tcrnbgh/Scratch/quarry_data/',
                 'prepData_alllocs_norm_maxdiff01_smpper50',
                 'int_times.RDS'))
print('done!')
# Clean up the cluster and release the relevant resources.
stopCluster(cl)
sink()
mpi.quit()