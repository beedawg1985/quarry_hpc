library(Rmpi)
library(snow)

setwd('/home/tcrnbgh/Scratch/quarry_data/quarry_hpc')
sink('./2_interpolate_analyse_hpc_sinkout.txt')

Sys.setenv(TMPDIR='/home/tcrnbgh/Scratch/tmp')
print('main R tempdir is...')
print(tempdir())

# prepare data --------------------------------
print('loading prepped data...')
f <- 'data/prepData_alllocs_norm_maxdiff01_smpper0.RDS'
prepData <- readRDS(f)
print('done!')
print('truncating prepData...')
prepData <- prepData[1:50]
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
  
  sessionTag <- 'prepData_alllocs_norm_maxdiff01_smpper0'
  
  sink(paste0('logs/2_interpolate_analyse_hpc_sinkout_site',pd$pol$fid,'.txt'))
  cvGrids <- loadCV()
  
  interpolateRas(pd,
    paramData = cvGrids,
    gLoc = grassGISDBASE,
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
save(datOut,paste0('/home/tcrnbgh/Scratch/quarry_data/',outputTag,'int_times.RDS'))
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
