#!/bin/bash -l
module purge
module load beta-modules
module load r/r-4.1.1_bc-3.13
module load nano/2.4.2
# required on compute nodes ?
module unload udunits
module load udunits/2.2.19
module load gerun

export PATH=/home/tcrnbgh/grass_latest/zstd/bin:$PATH
export LD_LIBRARY_PATH=/home/tcrnbgh/grass_latest/zstd/bin:$LD_LIBRARY_PATH

export PATH=/home/tcrnbgh/grass_latest/zstd/lib:$PATH
export LD_LIBRARY_PATH=/home/tcrnbgh/grass_latest/zstd/lib:$LD_LIBRARY_PATH

export GISBASE=/home/tcrnbgh/grass_latest/grass_install/grass80
 
export PATH=$PATH:$GISBASE/bin:$GISBASE/scripts
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GISBASE/lib

export PATH=$PATH:/home/tcrnbgh/grass_latest/grass_install/bin
export TZ="Europe/London"

