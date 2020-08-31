import pandas as pd
import matplotlib.cm as cm
import numpy as np
import matplotlib.pyplot as plt

def plot_clustered_stacked(dfall, labels=None, title="diBELLA V3 vs V1 (Overlap Only) —minus I/O",  H="/", **kwargs):
    """Given a list of dataframes, with identical columns and index, create a clustered stacked bar plot. 
labels is a list of the names of the dataframe, used for the legend
title is a string for the title of the plot
H is the hatch used for identification of the different dataframe"""

    n_df  = len(dfall)
    n_col = len(dfall[0].columns) 
    n_ind = len(dfall[0].index)
    axe = plt.subplot(111)
#  
    for df in dfall : # for each data frame
        axe = df.plot(kind="bar",
                      linewidth=0,
                      stacked=True,
                      ax=axe,
                      legend=False,
                      grid=False,
                      **kwargs)  # make bar plots

    h, l = axe.get_legend_handles_labels() # get the handles we want to modify
    for i in range(0, n_df * n_col, n_col): # len(h) = n_col * n_df
        for j, pa in enumerate(h[i:i+n_col]):
            for rect in pa.patches: # for each index
                rect.set_x(rect.get_x() + 1 / float(n_df + 1) * i / float(n_col))
                rect.set_hatch(H * int(i / n_col)) #edited part     
                rect.set_width(1 / float(n_df + 1))

    axe.set_xticks((np.arange(0, 2 * n_ind, 2) + 1 / float(n_df + 1)) / 2.)
    axe.set_xticklabels(df.index, rotation = 0)
    axe.set_title(title)

    # Add invisible data to add another legend
    n=[]        
    for i in range(n_df):
        n.append(axe.bar(0, 0, color="gray", hatch=H * i))

    l1 = axe.legend(h[:n_col], l[:n_col])

    if labels is not None:
        l2 = plt.legend(n, labels, loc=[1.01, 0.1]) 
    axe.add_artist(l1)
    return axe

# V3
df1io = pd.DataFrame(np.array([
                    [4.18,4.58,4.76,0.16,0.20,5.58,0.91,1.69,0.32+0.03],
                    [2.13,1.12,1.30,0.05,0.07,2.58,0.22,0.35,0.09+0.05],
                    [2.47,0.55,0.70,0.03,0.04,1.68,0.12,0.35,0.06+0.13],
                    [2.54,0.38,0.58,0.02,0.02,1.04,0.08,0.32,0.04+0.24]]),
                    index=["1", "4", "9", "16"],
                    columns=["DFD","KCFirstPass", "KCSecondPass", "SpMatA", "SpMatAt", "AAt", "ReadTaskExchange","WriteOverlap", "Other"])


# V1
df2io = pd.DataFrame(np.array([
                    [0.00,3.74,6.91,0.00,0.00,1.24,0.96,6.06,3.19+0.20],
                    [0.00,1.20,2.00,0.00,0.00,0.24,0.32,4.72,2.27+0.11],
                    [0.00,0.77,1.27,0.00,0.00,0.09,0.17,4.36,2.23+0.16],
                    [0.00,0.67,1.14,0.00,0.00,0.04,0.08,4.31,2.21+0.13]]),
                    index=["1", "4", "9", "16"],
                    columns=["DFD","KCFirstPass", "KCSecondPass", "SpMatA", "SpMatAt", "AAt", "ReadTaskExchange","WriteOverlap", "Other"])

# V3 no IO
df1 = pd.DataFrame(np.array([
                    [4.58,4.76,0.16,0.20,5.58,0.91,0.32+0.03],
                    [1.12,1.30,0.05,0.07,2.58,0.22,0.09+0.05],
                    [0.55,0.70,0.03,0.04,1.68,0.12,0.06+0.13],
                    [0.38,0.58,0.02,0.02,1.04,0.08,0.04+0.24]]),
                    index=["1", "4", "9", "16"],
                    columns=["KCFirstPass", "KCSecondPass", "SpMatA", "SpMatAt", "AAt", "ReadTaskExchange","Other"])

# V1 no IO
df2 = pd.DataFrame(np.array([
                    [3.74,6.91,0.00,0.00,1.24,0.96,3.19+0.20],
                    [1.20,2.00,0.00,0.00,0.24,0.32,2.27+0.11],
                    [0.77,1.27,0.00,0.00,0.09,0.17,2.23+0.16],
                    [0.67,1.14,0.00,0.00,0.04,0.08,2.21+0.13]]),
                    index=["1", "4", "9", "16"],
                    columns=["KCFirstPass", "KCSecondPass", "SpMatA", "SpMatAt", "AAt", "ReadTaskExchange","Other"])

#################### intra node on 4 node

# V3
df1io4 = pd.DataFrame(np.array([
                    [1.03,70.38,77.49,2.62,2.91,28.95,40.02,41.43, 4.60+0.02],
                    [0.89,18.02,18.02,0.67,0.71,12.40,5.44,5.79,   1.16+0.07],
                    [0.91,4.32,4.16,0.14,0.19,5.67,0.89,1.17,      0.31+0.05],
                    [2.13,1.12,1.30,0.05,0.07,2.58,0.22,0.35,      0.09+0.05]]),
                    index=["1", "4", "16", "64"],
                    columns=["DFD","KCFirstPass", "KCSecondPass", "SpMatA", "SpMatAt", "AAt", "ReadTaskExchange","WriteOverlap", "Other"])

# V1
df2io4 = pd.DataFrame(np.array([
                    [0.00,53.17,81.28,0.00,0.00,30.54,13.56,8.68,6.80+0.38],
                    [0.00,13.21,19.23,0.00,0.00,6.18,3.35,6.29,  6.11+0.17],
                    [0.00,3.39,7.25,0.00,0.00,1.21,0.89,4.73,2.34+0.23],
                    [0.00,1.20,2.00,0.00,0.00,0.24,0.32,4.72,2.27+0.11]]),
                    index=["1", "4", "16", "64"],
                    columns=["DFD","KCFirstPass", "KCSecondPass", "SpMatA", "SpMatAt", "AAt", "ReadTaskExchange","WriteOverlap", "Other"])

# V3 no io
df14 = pd.DataFrame(np.array([
                    [70.38,77.49,2.62,2.91,28.95,40.02, 4.60+0.02],
                    [18.02,18.02,0.67,0.71,12.40,5.44, 1.16+0.07],
                    [4.32,4.16,0.14,0.19,5.67,0.89,  0.31+0.05],
                    [1.12,1.30,0.05,0.07,2.58,0.22, 0.09+0.05]]),
                    index=["1", "4", "16", "64"],
                    columns=["KCFirstPass", "KCSecondPass", "SpMatA", "SpMatAt", "AAt", "ReadTaskExchange","Other"])

# V1 no io
df24 = pd.DataFrame(np.array([
                    [53.17,81.28,0.00,0.00,30.54,13.56,6.80+0.38],
                    [13.21,19.23,0.00,0.00,6.18,3.35, 6.11+0.17],
                    [3.39,7.25,0.00,0.00,1.21,0.89,2.34+0.23],
                    [1.20,2.00,0.00,0.00,0.24,0.32,2.27+0.11]]),
                    index=["1", "4", "16", "64"],
                    columns=["KCFirstPass", "KCSecondPass", "SpMatA", "SpMatAt", "AAt", "ReadTaskExchange", "Other"])

# Then, just call 
plt.figure(figsize=(11, 6))
plt.rcParams.update({'font.size': 14})

# plt.ylabel('Time (s)')
# plt.xlabel('Processors (64 MPI Rank/Processor')
plt.ylabel('Time (s)')
plt.xlabel('Processors (64 MPI Rank/Processor)')

plot_clustered_stacked([df1,df2],["V3", "V1"], cmap=plt.cm.tab20)

plt.show()