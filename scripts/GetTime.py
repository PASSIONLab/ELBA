# https://pythonicways.wordpress.com/2016/12/20/log-file-parsing-in-python/
# Export parsed data to text file

import re
import time
import sys
from time import strftime
 
LogFileName = open(sys.argv[1], 'r')

# @GGGG: run the program from build folder (temporary)
ExportFile = "GetTime-" + str(sys.argv[1]) + ".txt"

regex1st = '(Main|\(\):[0-9]).*' # Only care about the last lines of the log
regex2nd = '[0-9].*[0-9]'        # Extract only numeric values 
ReadLine = True
ToSecond = True
 
with open(sys.argv[1], "r") as file:
    MatchList = []
    for line in file:
        for match in re.finditer(regex1st, line, re.S):
            MatchText = match.group()
            for num in re.finditer(regex2nd, MatchText, re.S):
                MatchText2 = num.group()
                MatchList.append(MatchText2)
                # print MatchText2
file.close()
 
with open(ExportFile, "w+") as file:
    for item in [x for x in xrange(0, len(MatchList)) if x != 2 and x != 3 and x != 8]: # Remove lines I'm not interested in
        if ToSecond == True:
            ms = float(MatchList[item]) # Trasform from ms to s
            s  = ms / 1000
            print str(s)
            file.write(str(s) + "\n")
        else:
            print MatchList[item]
            file.write(MatchList[item] + "\n")
file.close()
