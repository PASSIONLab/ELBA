# https://stackoverflow.com/questions/208120/how-to-read-and-write-multiple-files
# Read multiples files and compute avg time among them

import re
import time
import sys
import glob
import numpy as np
 
LogFileName = str(sys.argv[1])
ExportFile  = "/global/cscratch1/sd/gguidi/IPDPS2021/performance/"
ExportFile  = ExportFile + "ComputeStats-" + LogFileName + ".txt" 

sums = np.zeros(12)

# @GGGG: run this within a folder with the files I'm interested in averaging
ListOfFiles = glob.glob('GetTime-' + LogFileName + '*.txt') # Create the list of file that I already parsed
for ThisFile in ListOfFiles:
  FI = open(ThisFile, 'r')
  numbersFromFile = [float(line) for line in FI.readlines()]
  array = np.array(numbersFromFile)
  sums  = sums + array

# Compute avg times
sums = sums / len(ListOfFiles)

# Print results and write them to file
with open(ExportFile, "w+") as file:
    file.write(LogFileName + "\n")  # First line is gonna be a header
    for i in xrange(0, len(sums)):
        print sums[i]
        file.write(str(sums[i]) + "\n")
file.close()
