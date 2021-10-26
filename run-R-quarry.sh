#!/bin/bash -l

#  R MPI parallel job

# Request ten minutes of wallclock time (format hours:minutes:seconds).
#$ -l h_rt=02:00:0

# Request 1 gigabyte of RAM per process.
#$ -l mem=2G

# Set the name of the job.
#$ -N quarry_interpolations_test

# Select the MPI parallel environment with 32 processes
#$ -pe mpi 70

# Set the working directory to somewhere in your scratch space.  This is
# necessary because the compute nodes cannot write to your $HOME
# NOTE: this directory must exist.

# Load R/GRASS environment
source /home/tcrnbgh/quarry_hpc/init.sh

# run script to make working directory
source /home/tcrnbgh/quarry_hpc/mkwd.sh

# set working dir
#$ -wd /home/tcrnbgh/Scratch/quarry_data

# Run our MPI job. GERun is our wrapper for mpirun, which launches MPI jobs  
gerun RMPISNOW < /home/tcrnbgh/Scratch/quarry_data/quarry_hpc/rscript/2_interpolate_analyse_hpc.R > 
snow.out.${JOB_ID}


