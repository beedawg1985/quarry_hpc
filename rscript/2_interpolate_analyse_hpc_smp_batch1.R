library(Rmpi)
library(snow)

setwd('/home/tcrnbgh/Scratch/quarry_data/quarry_hpc')
sink('./2_interpolate_analyse_hpc_sinkout.txt')

Sys.setenv(TMPDIR='/home/tcrnbgh/Scratch/tmp')
print('main R tempdir is...')
print(tempdir())

# prepare data --------------------------------
print('loading prepped data...')
f <- 'data/prepData_alllocs_norm_maxdiff01_hole2.RDS'             # check this !
prepData <- readRDS(f)
length(prepData)

print('done!')
print('truncating prepData...')
prepData <- prepData[1:160]                                       # check this !

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

  sessionTag <- 'prepData_alllocs_norm_maxdiff01_hole2'           # check this !
  
  sink(paste0('logs/2_interpolate_analyse_hpc_sinkout_site',
              pd$pol$fid,'.txt'))
  
  cvGrids <- loadCV2()                                            # check this !
  
  interpolateRas(pd,
    paramData = cvGrids,
    gLoc = grassGISDBASE,
    outputDir = '/home/tcrnbgh/Scratch/quarry_data/data_output',
    testCV = F, # = T for test run                                # check this !
    outputTag = sessionTag,
    intMethods=c(
      'rfsp',
      'nn','idw','ok','tin',
      'gbicubic',
      'gspline'
      )
    )
  
})
save(datOut,
     file=paste0('/home/tcrnbgh/Scratch/quarry_data/',
                 'prepData_alllocs_norm_maxdiff01_hole2',
                 'int_times.RDS'))
print('done!')
# Clean up the cluster and release the relevant resources.
stopCluster(cl)
sink()
mpi.quit()