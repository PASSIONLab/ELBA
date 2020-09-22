# https://ask.sagemath.org/question/42866/how-can-i-generate-a-graph-from-an-mtx-file/

import sys
import scipy
from scipy.io import mmread
# from graphviz import Digraph
import math 

# For SORA inpute, an edge has 9 attributes: 3,F,33,0,0,2,34,0,32
# Col1: overlap orientation [@GGGG: check if these are consistent with the ones I defined]
# 0 = u<--------<v      reverse of u to reverse of v    [@GGGG: equivalent to 2 <--< in diBELLA]
#   => This case is handled in DOT file preprocessing step and changed to 3 (u>-->v)
# 1 = u<-------->v      reverse of u to forward of v    [@GGGG: equivalent to 3 <--> in diBELLA]
# 2 = u>--------<v      forward of u to reverse of v    [@GGGG: equivalent to 0 >--< in diBELLA]
# 3 = u>-------->v      forward of u to forware of v    [@GGGG: equivalent to 1 >--> in diBELLA]
# Col2: overlap property F:forward, [@GGGG: rc = 0] 
#                        FRC:read1 overlaps with the reverse complement of read2 [@GGGG: rc = 1]
# Col3~9: overlap length, substitutions, edits, start1, stop1, start2, stop2 [@GGGG: substitutions, edits are gonna be always zero]
# digraph G {
# 	7.1
# 	7.1 -> 2.1 [label="F,31,0,0,35,0,30,35,4,34"]
# 	7.1 -> 3.1 [label="F,34,0,0,35,0,33,35,1,34"]
# 	8.1
# 	8.1 -> 7.1 [label="F,34,0,0,35,0,33,35,1,34"]
# 	8.1 -> 2.1 [label="F,30,0,0,35,0,29,35,5,34"]
#   ...
# 	6.1
# 	6.1 -> 4.1 [label="F,33,0,0,35,0,32,35,2,34"]
# 	6.1 -> 3.1 [label="F,32,0,0,35,0,31,35,3,34"]
# }

# This is diBELLA's AAt output format [@GGGG: need to output direction, rc, ovlen, start and end positions]
# direction, rc, begV, endV, begH, endH (OverlapLen and others computed in python script during translation)
# %MatrixMarket matrix coordinate real general
# 16890 16890 756
# ...

# Read diBELLA input
InputFileName  = open(sys.argv[1], 'r')
ReadFileName   = open(sys.argv[2], 'r')

OutputFileName = str(sys.argv[1])

lines = InputFileName.readlines()
names = ReadFileName.readlines()

del lines[0] # Remove mmx header
del lines[0] # Remove line with summary nodes/edges

# @GGGG: create edge list format (SORA) and paf format (miniasm)

# SORA
# 0	4	3,F,32,0,0,3,34,0,31
# 2	4	3,F,33,0,0,2,34,0,32
# 0	2	3,F,34,0,0,1,34,0,33

edgeslist = open(OutputFileName + "_edge_list.txt", "w")
for line in lines:
    items = line.split("\t")

    vertex1 = items[0]
    vertex2 = items[1]
    
    direction = items[2]

    if(direction == '0'): 
        direction = '2'         # <-->
    elif(direction == '1'):
        direction = '3'         # >-->
    elif(direction == '2'):
        direction = '0'         # <--<
    elif(direction == '3'):
        direction = '1'         # <-->

    rc = items[3]

    if(rc == '1'):
        rc = "FRC"
    else:
        rc = "F"

    begV = items[5]
    endV = items[6]
    begH = items[7]
    endH = items[8]

    overlap = items[11]

    entry = str(vertex1) + '\t' + str(vertex2) + '\t' + str(direction) + ',' + rc + ',' + str(overlap) + ',' + str(0) + ',' + str(0) + ',' + str(begV) + ',' + str(endV) + ',' + str(begH) + ',' + str(endH) + '\n'
    edgeslist.write(entry)

edgeslist.close()

# miniasm 
# 1	string	Query sequence name
# 2	int	Query sequence length
# 3	int	Query start (0-based; BED-like; closed)
# 4	int	Query end (0-based; BED-like; open)
# 5	char	Relative strand: "+" or "-"
# 6	string	Target sequence name
# 7	int	Target sequence length
# 8	int	Target start on original strand (0-based)
# 9	int	Target end on original strand (0-based)
# 10	int	Number of residue matches
# 11	int	Alignment block length
# 12	int	Mapping quality (0-255; 255 for missing)

def toOriginalCoordinates(begpH, endpH, lenH):
    tmp   = begpH
    begpH = lenH-endpH
    endpH = lenH-tmp
    return begpH, endpH

# first, parse readNameMap and create a dictionary
readnamemap = {}
for name in names:
    items = name.split("\t")
    idx   = items[0]
    name  = items[1]
    name  = name[1:]
    readnamemap[int(idx)] = name 

def ReadName(idx):
    return readnamemap[idx]

pafformat = open(OutputFileName + "_miniasm.paf",   "w")
for line in lines:
    items = line.split("\t")

    vertex1 = items[0]
    vertex2 = items[1]

    nameV = ReadName(int(vertex1))
    nameH = ReadName(int(vertex2))

    nameV = nameV.rstrip("\n")
    nameH = nameH.rstrip("\n")

    rc = items[3]

    begV = items[5]
    endV = items[6]
    begH = items[7]
    endH = items[8]

    lenV = items[9]
    lenH = items[10]

    if(rc == '1'):
        rc = "-"
        begH, endH = toOriginalCoordinates(int(begH), int(endH), int(lenH))
    else:
        rc = "+"

    alnlen     = items[11]
    resmatches = math.floor(0.8*float(alnlen))
    mapping    = '0'

    entry = nameV + '\t' + str(lenV) + '\t' + str(begV) + '\t' + str(endV) + '\t' + rc + '\t' + nameH + '\t' + str(lenH) + '\t' + str(begH) + '\t' + str(endH) + '\t' + str(resmatches) + '\t' + alnlen + '\t' + mapping + '\n'
    pafformat.write(entry)

pafformat.close()

ReadFileName.close()
InputFileName.close()

# Create a digraph with possible parallel edges and self-loops use
# dot = Digraph('G')

# # @GGGG: parse the input file and populate the graph 
# for line in lines:
#     items = line.split("\t")

#     vertex1 = items[0]
#     vertex2 = items[1]
    
#     direction = items[2]
#     rc        = items[3]

#     if(rc == '1'):
#         rc = "FRC"
#     else:
#         rc = "F"

#     begV = items[4]
#     endV = items[5]
#     begH = items[6]
#     endH = items[7]
#     overlap   = 0 #items[8]

#     currlabel = str(direction) + ',' + rc + ',' + str(overlap) + ',' + str(0) + ',' + str(0) + ',' + str(begV) + ',' + str(endV) + ',' + str(begH) + ',' + str(endH)

#     dot.node(vertex1)  # Create a node only for the first read
#     dot.edge(vertex1, vertex2, label=currlabel)  # Create an edge between the first and second node and populate the label with the overlap information