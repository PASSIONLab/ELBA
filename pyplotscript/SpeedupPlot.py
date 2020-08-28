import matplotlib
import matplotlib.pyplot as plt
import numpy as np

v3_overall = (2.8,3.7,4.2)
v1_overall = (2.1,2.5,2.6)

v3_overlap = (2.2,3.3,5.4)
v1_overlap = (5.2,13.8,31.0)

v3_exchange = (4.1,7.5,10.7)
v1_exchange = (3.0,5.4,12.1)

ind = np.arange(len(v1_overall))  # the x locations for the groups
width = 0.15  # the width of the bars

plt.rcParams.update({'font.size': 14})

fig, ax = plt.subplots(figsize=(11, 6))

rects1 = ax.bar(ind - 0.375, v3_overall, width,
                label='V3 Overall', color='#88498F')
rects2 = ax.bar(ind - 0.225, v1_overall, width,
                label='V1 Overall', color='#779FA1')
rects3 = ax.bar(ind - 0.075, v3_overlap, width,
                label='V3 AAt',color='#E0CBA8')
rects4 = ax.bar(ind + 0.075, v1_overlap, width,
                label='V1 AAt',color='#FF6542')
rects5 = ax.bar(ind + 0.225, v3_exchange, width,
                label='V3 ReadTaskExchange',color='#AFC97E')
rects6 = ax.bar(ind + 0.375, v1_exchange, width,
                label='V1 ReadTaskExchange',color='#E2D686')

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Speed-Up')
ax.set_xlabel('Processors (64 MPI Rank/Processor)')
ax.set_title('diBELLA V3 vs V1 Speed-Up (Overlap Only)')
ax.set_xticks(ind)
ax.set_xticklabels(("4", "9", "16"))
ax.legend(loc=2)

ax.axhline(y=4, linewidth=1, color='r', linestyle='--')
ax.axhline(y=9, linewidth=1, color='r', linestyle='-.')
ax.axhline(y=16, linewidth=1, color='r',linestyle=':')

def autolabel(rects, xpos='center'):
    """
    Attach a text label above each bar in *rects*, displaying its height.

    *xpos* indicates which side to place the text w.r.t. the center of
    the bar. It can be one of the following {'center', 'right', 'left'}.
    """

    ha = {'center': 'center', 'right': 'left', 'left': 'right'}
    offset = {'center': 0, 'right': 1, 'left': -1}

    for rect in rects:
        height = rect.get_height()
        ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(offset[xpos]*3, 3),  # use 3 points offset
                    textcoords="offset points",  # in both directions
                    ha=ha[xpos], va='bottom')

# autolabel(rects1, "left")
# autolabel(rects2, "right")

fig.tight_layout()

plt.show()