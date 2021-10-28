#!/usr/bin/env python
import csv
import sys
import os
import datetime
import grass.script as grass
import multiprocessing as multi
import itertools

def my_print(text):
    currentDT = datetime.datetime.now()
    sys.stdout.write('\n')
    sys.stdout.write(str(currentDT))
    sys.stdout.write('\n' + str(text))
    sys.stdout.flush()

def runInt(smoothVal,tensionVal,npminVal,runNo,polfid,outputDir):
  my_print('running v.surf.rst...')
  grass.run_command("v.surf.rst", input='points',\
                    segmax='10',\
                    zcolumn='elev',\
                    tension=smoothVal,\
                    smooth=tensionVal,\
                    npmin=npminVal,\
                    nprocs='1',\
                    mask='testMask',\
                    elevation='int',\
                    overwrite=True)
  grass.run_command("r.out.gdal",\
                    input='int',\
                    output=outputDir + '/gspline_int_intfid_' + \
                    polfid + '_runnum_' + runNo + '.tif',
                    overwrite=True)
  return
   
# parse arguments
smoothVal = sys.argv[1]
tensionVal = sys.argv[2]
npminVal = sys.argv[3]
runNo = sys.argv[4]
polfid = sys.argv[5]
outputDir = sys.argv[6]
# runInt(2,1e-04,150,2233,68,'/lustre/scratch/tcrnbgh/quarry_data/quarry_hpc/raster')
attempts = 0
while attempts < 6:
  try:
    runInt(smoothVal,tensionVal,npminVal,runNo,polfid,outputDir)
    break
  except:
    my_print('error running interpolation, trying again...')
    attempts += 1
