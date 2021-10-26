#!/usr/bin/env Rscript
setwd('/home/tcrnbgh/Scratch/quarry_data/quarry_hpc')

# !!!!
# have you run cp ./quarry_hpc/rscript/* ./Scratch/quarry_data/quarry_hpc/rscript/
# system('cp /home/tcrnbgh/quarry_hpc/rscript/* /home/tcrnbgh/Scratch/quarry_data/quarry_hpc/rscript/')
# to copy repo to scratch directory?

# check if packages need installing
list.of.packages <- c("tuneRanger","gstat","reshape2",
                      "e1071","caret","randomForest","stringr",
                      "stars","sf","dplyr","gdalUtils",
                      "raster","automap","fields","interp",
                      "mgcv","purrr","furrr","doParallel",
                      "future.apply", "snow")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages,
                                          lib="/home/tcrnbgh/R/x86_64-pc-linux-gnu-library/4.1")

# set envs
userDataDir <<- '/home/tcrnbgh/Scratch/quarry_data'
grassMapset <<- paste0(userDataDir,'/grassdb/quarry/PERMANENT/')
grassLocation <<- paste0(userDataDir,'/grassdb/quarry/')

# source functions
source('rscript/general_functions.R')

# new CV grids ----------------------------------------------------
# prep parameters for grass resamp.filter 
f1 = c( 'gauss', 'normal', 'sinc', 'hann', 'hamming', 'blackman')
f2 = c('hermite',  'lanczos1', 'lanczos2', 'lanczos3','box', 'bartlett')
df <- rbind(data.frame(expand.grid(f1,f2)),data.frame(Var1 = f2, Var2 = NA))
gfilter.params <- df %>% 
  tidyr::unite('col3', c(Var1,Var2), sep = ',', na.rm = TRUE) %>% 
  mutate_if(is.factor,as.character) %>% 
  mutate(filtcomb_number = 1:nrow(.))
save(gfilter.params,file='cvdev/gfilter.RDS')

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

# end new CV grids ------------------------------------------------

# or use existing prepared data --------------------------------

f <- 'data/prepData_alllocs_maxdiff01_smpper0.RDS'
prepData <- readRDS(f)

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
cl <- getMPIcluster()
#cl <- snow::makeMPIcluster(count=70)
# cl <- makeCluster(mc <- getOption("cl.cores", 60), outfile='logs/cluster_out.txt')
clusterEvalQ(cl, source(paste0(getwd(),'/rscript/general_functions.R')))
snow::clusterExport(cl, list=c('cvGrids','grassLocation','sessionTag',
                               'userDataDir','grassMapset'))

datOut <- parLapply(cl, prepData[1:70], function(pd) {
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
