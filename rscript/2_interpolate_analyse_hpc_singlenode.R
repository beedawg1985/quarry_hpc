library(snow)

setwd('/home/tcrnbgh/Scratch/quarry_data/quarry_hpc')
sink('./2_interpolate_analyse_hpc_sinkout.txt')

Sys.setenv(TMPDIR='/home/tcrnbgh/Scratch/tmp')
print('main R tempdir is...')
print(tempdir())

# prepared data --------------------------------
print('loading prepped data...')
f <- "data/prepData_alllocs_norm_maxdiff01_smpper50.RDS"
prepData <- readRDS(f)
print('done!')
print('truncating prepData...')
prepData <- prepData[1:35]
print('done!')
# interpolation run --------------------------------------------
print('running interpolations...')
st <- Sys.time()

#!! if offSet = T the use pd$tiles$pol !!
print('making node cluster...')
cl <- parallel::makeCluster(10,
                            type='PSOCK')
# print('done!')
print('running interpolations...')

# pd <- prepData[[1]]
datOut <- snow::clusterApplyLB(cl, prepData, function(pd) {
  
  setwd('/home/tcrnbgh/Scratch/quarry_data/quarry_hpc')
  source('rscript/interpolation_functions.R')
  print('temp dir...')
  print(tempdir())
  Sys.setenv(TMPDIR='/home/tcrnbgh/Scratch/tmp')
  
  userDataDir <- '/home/tcrnbgh/Scratch/quarry_data'
  grassGISDBASE <- paste0(userDataDir,'/grassdb')
  
  sessionTag <- 'prepData_alllocs_norm_maxdiff01_smpper0'
  
  sink(paste0('logs/2_interpolate_analyse_hpc_sinkout_site',
              pd$pol$fid,'.txt'))
  
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
                   'gbicubic',
                   'gspline'
                 )
  )
  
})
save(datOut,file=paste0('/home/tcrnbgh/Scratch/quarry_data/',outputTag,'int_times.RDS'))
print('done!')
# Clean up the cluster and release the relevant resources.
stopCluster(cl)
sink()
mpi.quit()

