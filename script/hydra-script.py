#!/usr/bin/python
# -*- coding: utf-8 -*-

# G. Guidi -- Compute network degree and neighborhood size (ELBA)

from scipy.io import mmread
from scipy.sparse import csr_matrix
import numpy as np
import sys
import networkx as nx
import pandas as pd
import matplotlib.cm as cm
import matplotlib.font_manager
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

def readspmat(filename, dtype: str = "int32"):
    """
    Read `.mm` file.

    Parameters
    ----------
    filename
        The filename.
    dtype
        Numpy data type.
    """
    X = mmread(filename).astype(dtype)
    X = csr_matrix(X)

    return X

def degreedistribution(m):
    
    # nnz in my matrix
    nnz = m.nnz
    print("nnz: ", nnz)

    G = nx.from_scipy_sparse_matrix(m) # iterate over the matrix
    nc = G.number_of_nodes() # node count

    # the array index is the size of the neighborhood and the array value is the number of sequences having that neighborhood size
    d1n = np.zeros(nc)
    d2n = np.zeros(nc)
    d3n = np.zeros(nc)
    d4n = np.zeros(nc)

    dc = 0 # degree count -- shold be == nnz?
    hd = 0 # high degree vertex

    maxd1 = 0 # max 1-ring neighborhood size
    maxd2 = 0 # max 1-ring neighborhood size
    maxd3 = 0 # max 1-ring neighborhood size
    maxd4 = 0 # max 1-ring neighborhood size

    for node in range(nc):
        degree = G.degree[node]
        d1n[degree] += 1

        dc += degree

        if degree > maxd1:
            maxd1 = degree

    print("max 1-ring neighborhood size: ", maxd1)
    print("max 2-ring neighborhood size: ", maxd2)
    print("max 3-ring neighborhood size: ", maxd3)
    print("max 4-ring neighborhood size: ", maxd4)

    plot(x, y, 'bo')  # plot x and y using blue circle markers

    return d1n, d2n, d3n, d4n

def main():
    m = readspmat(sys.argv[1])
    d1, d2, d2, d4 = degreedistribution(m)
    print(d1)
    print("End of program ;)")

if __name__ == "__main__":
    main()