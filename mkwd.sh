#!/bin/bash -l
# remove and remake working dir
rm -rf /home/tcrnbgh/Scratch/quarry_data 
mkdir /home/tcrnbgh/Scratch/quarry_data

# run grass to initialise location
grass -c /home/tcrnbgh/quarry_hpc/vector/intout_1_training.gpkg /home/tcrnbgh/Scratch/quarry_data/grassdb/quarry -e

export GISDBASE=/home/tcrnbgh/Scratch/quarry_data/grassdb

# create output folder
mkdir /home/tcrnbgh/Scratch/quarry_data/data_output

# copy R repo
cp -r /home/tcrnbgh/quarry_hpc /home/tcrnbgh/Scratch/quarry_data/quarry_hpc

# make python scripts executable
chmod +x /home/tcrnbgh/Scratch/quarry_data/quarry_hpc/python/GRASS_vrst.py
chmod +x /home/tcrnbgh/Scratch/quarry_data/quarry_hpc/python/GRASS_bspline.py
chmod +x /home/tcrnbgh/Scratch/quarry_data/quarry_hpc/python/GRASSoptimise.py
chmod +x /home/tcrnbgh/Scratch/quarry_data/quarry_hpc/python/GRASS_resampfilter.py