# library(Rmpi)
library(snow)
library(snowfall)
library(stringr)

setwd('/home/tcrnbgh/Scratch/quarry_data/quarry_hpc')
sink('./2_interpolate_analyse_hpc_sinkout.txt')

Sys.setenv(TMPDIR='/home/tcrnbgh/Scratch/tmp')
print('main R tempdir is...')
print(tempdir())

# prepare data --------------------------------
print('loading prepped data...')
f <- 'data/prepData_alllocs_norm_maxdiff01_hole2.RDS'             # check this !
prepData <- readRDS(f)
print('done!')

print('truncating prepData...')
prepData <- prepData[1:60]                                       # check this !
print('done!')

# interpolation run --------------------------------------------
print('running interpolations...')
st <- Sys.time()

#!! if offSet = T the use pd$tiles$pol !!

print('processing tasks...')
methods <- list('rfsp','nn','ok','tin',
                'idw','gbicubic','gspline')
tasks <- expand.grid(methods,1:60)                                # check this !
sessionTag <- str_replace(basename(f),'.RDS','')
print('done!')


print('making node cluster...')
snowfall::sfInit(parallel=T, cpus=60, type='MPI')                 # check this !
snowfall::sfInit(parallel=T, cpus=30, type='MPI')                 # check this !
snowfall::sfExport(list=c('prepData','tasks','sessionTag'))
print('done!')

print('running interpolations...')

datOut <- snowfall::sfClusterApplyLB(1:nrow(tasks), function(tnum) {
  setwd('/home/tcrnbgh/Scratch/quarry_data/quarry_hpc')
  source('rscript/interpolation_functions.R')
  print('temp dir...')
  print(tempdir())
  Sys.setenv(TMPDIR='/home/tcrnbgh/Scratch/tmp')
  
  userDataDir <- '/home/tcrnbgh/Scratch/quarry_data'
  grassGISDBASE <- paste0(userDataDir,'/grassdb')
  
  sink(paste0('logs/2_interpolate_analyse_hpc_hole_lb_sinkout_site',
              Sys.getpid(),'.txt'))
  
  m <- tasks[tnum,'Var1']
  interpolateRas(prepData[[tasks[tnum,'Var2']]],
                  paramData=loadCV(),
                  gLoc = 
                    '/home/tcrnbgh/Scratch/quarry_data/grassdb',
                  outputDir = 
                    '/home/tcrnbgh/Scratch/quarry_data/data_output',
                  testCV = T,                                     # check this !
                  outputTag = 
                    paste0(sessionTag,'_',paste0(m,collapse='_')),
                  intMethods=m
  )
})


save(datOut,
     file=paste0('/home/tcrnbgh/Scratch/quarry_data/',
                 'prepData_alllocs_norm_maxdiff01_hole2',
                 'int_times.RDS'))
print('done!')

# Clean up the cluster and release the relevant resources.
sink()
mpi.quit()