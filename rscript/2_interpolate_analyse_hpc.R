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

.GlobalEnv$.MPIrun = MPIrun()
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
if (is.null(cl)) { 
  print('getMPIcluster() failed...')
  print('trying to make cluster...')
  cl <- 
    snow::makeMPIcluster(count=35,
                         outfile=paste0(getwd(),'/logs/cluster_out.txt'))
  cl <- makeCluster(mc <- getOption("cl.cores", 35),
                    outfile=paste0(getwd(),'/logs/cluster_out.txt'))
  }
print('done!')

print('sourcing functions and exporting vars on cluster...')
clusterEvalQ(cl, source(paste0(getwd(),'/rscript/general_functions.R')))
snow::clusterExport(cl, list=c('cvGrids','grassLocation','sessionTag',
                               'userDataDir','grassMapset'))
print('done!')

print('running interpolations...')
prepDataTrunc <- prepData[1:70]
datOut <- parLapply(cl, prepDataTrunc, function(pd) {
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
