#!/bin/bash -l

#  R MPI parallel job

# Request ten minutes of wallclock time (format hours:minutes:seconds).
#$ -l h_rt=13:00:0

# Request 1 gigabyte of RAM per process.
#$ -l mem=300M

# Set the name of the job.
#$ -N quarry_interpolations_batch3

# Select the MPI parallel environment with 32 processes
#$ -pe mpi 160

# Load R/GRASS environment
echo "running init.sh script..."
source /home/tcrnbgh/quarry_hpc/init.sh
echo "done!"

# Set the working directory to somewhere in your scratch space.  This is
# necessary because the compute nodes cannot write to your $HOME
# NOTE: this directory must exist.
# set working dir
#$ -wd /home/tcrnbgh/Scratch/quarry_data

# Run our MPI job. GERun is our wrapper for mpirun, which launches MPI jobs  
echo "running gerun..."
gerun /home/tcrnbgh/quarry_hpc/RMPISNOW_bgh < /home/tcrnbgh/Scratch/quarry_data/quarry_hpc/rscript/2_interpolate_analyse_hpc_batch3.R
echo "done!"

