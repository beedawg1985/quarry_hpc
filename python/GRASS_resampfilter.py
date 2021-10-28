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

def runInt(radVal,filterVal,runNo,polfid,outputDir):
  my_print('running r.resamp.filter...')
  grass.run_command("r.resamp.filter", input='allRas',\
                    output='filterRas',\
                    filter_=filterVal,\
                    radius=radVal,\
                    overwrite=True)
  grass.run_command("r.out.gdal",\
                    input='filterRas',\
                    output=outputDir + '/gfilter_int_intfid_' + \
                    polfid + '_runnum_' + runNo + '.tif',
                    overwrite=True)
  return
   
# parse arguments
radVal = sys.argv[1]
filterVal = sys.argv[2]
runNo = sys.argv[3]
polfid = sys.argv[4]
outputDir = sys.argv[5]

while attempts < 6:
  try:
    runInt(radVal,filterVal,runNo,polfid,outputDir)
    break
  except:
    my_print('error running interpolation, trying again...')
    attempts += 1


