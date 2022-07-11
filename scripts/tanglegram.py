import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection, LineCollection
import matplotlib.path as mpath
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

import numpy as np
import pandas as pd

import itertools
import time

import xio

import sys

sys.path.append("./baltic3-master")

import baltic3 as bt
import baltic3_utils as btu


def get_leaf_coords(tipname, tre):
    """Searches for part of, or all of, a tipname in a tree. 
    Returns the first instance if multiple matches exist, but this fails silently. """

    found = 0
    tip_x = 0
    tip_y = 0
    for lf in tre.leaves:
        if tipname in lf.name:
            tip_x = lf.height
            tip_y = lf.y
            found = 1
    if found == 0:
        print("ERROR: %s not found!" % tipname)
        return -1, -1
    return tip_x, tip_y


#import data
seg_ls = ["HA","NA","PB2", "PB1", "PA", "NP", "MP", "NS"]

df_dict = {}
for seg in seg_ls:
    df_dict[seg] = btu.austechia_read_tree("../data/Tanglegram/"+seg+"_sub.afa.nex")

df = pd.read_csv('../data/Tanglegram/HA_sub_clade_assignment.tsv', sep='\t', header=0)
df['Taxa'] = df['Taxa'].replace("\\|.*","",regex=True)
df = df.set_index('Taxa').T.to_dict('records')[0]
    
    
# Plot tanglegram!
# =================================== PARAMS ===================================

#zoom
y_zoom_dict = {}
y_zoom_dict["HA"] = 1.45
y_zoom_dict["NA"] = 1.27
y_zoom_dict["NS"] = 1.06
y_zoom_dict["NP"] = 1.05
y_zoom_dict["MP"] = 1.145
y_zoom_dict["PA"] = 1.015
y_zoom_dict["PB1"] = 1.005
y_zoom_dict["PB2"] = 1
# for seg in seg_ls:
#     #if seg == "HA":
#     #    y_zoom_dict[seg] = 1.35
#     #else:
# 	y_zoom_dict[seg] = 1

# TREE PARAMS
branchWidth=0.5 # line thickness of branches
dotted_line_width=0.5
dotted_line_alpha=0.3
global_ySpan = df_dict["HA"].ySpan/y_zoom_dict["HA"]

# x coordinates to offset each tree by
# key: val = seg: [cumulative+tree.Ht, cumulative+space]
x_offset_dict = {}
cumulative_x_offset = 0
for seg in seg_ls:
    print(seg)
    print(cumulative_x_offset)
    x_offset_dict[seg] = cumulative_x_offset
    cumulative_x_offset += df_dict[seg].treeHeight*1.5


y_offset_dict = {}
for seg in seg_ls:  
    y_offset_dict[seg] = 0



# Prep next seg dict
next_seg_dict = {}

for i in range(len(seg_ls)-1):
    next_seg_dict[seg_ls[i]] = seg_ls[i+1]


t0 = time.time()

fig,ax = plt.subplots(figsize=(20, 9.3),facecolor='w')
patches_ls = []
lines_ls = []


cdict = {"3a1": "#40ab5d", 
         '3a2': "#e06c9a", 
         '3a2 (global)': "#fed439",
         "3a": "#e97c71",
         "V1A.3": "#7f7dbb", 
         "V1A": "#6badd6",
         "unassigned": "black"
        }
  
for i in range(len(seg_ls)):
    seg = seg_ls[i]
    
    x_offset = x_offset_dict[seg]
    y_offset = y_offset_dict[seg]
    print(seg)
    

    ax.text(x_offset+df_dict[seg].treeHeight*1.5/4, -30, seg, size=20)
    
    for k in df_dict[seg].Objects:
        c = 'k'
        x=k.height
        y=k.y/y_zoom_dict[seg]

        xp = k.parent.height
        if x is None: # matplotlib won't plot Nones, like root
            x = x_offset
        if xp==None:
            xp = x + x_offset

        if isinstance(k,bt.leaf) or k.branchType=='leaf':
            if seg != seg_ls[-1]:
                # Coords of the current tip
                x0 = x+x_offset
                y0 = y+y_offset
                
                # Get coords of the tip in the next tree
                if seg == "HA":
                	tip_ident = k.name.split("|")[0]
                else:
                	tip_ident = k.name
                print(tip_ident)
                clade_label = df[tip_ident]
                if clade_label == "3a2" or clade_label == "3a2 (global)" or clade_label == "V1A":
                    dotted_line_alpha = 0.7
                else:
                    dotted_line_alpha = 0.2
				           
                next_x_offset = x_offset_dict[next_seg_dict[seg]]
                next_y_offset = y_offset_dict[next_seg_dict[seg]]
                x1, y1 = get_leaf_coords(tip_ident, df_dict[next_seg_dict[seg]])
                
                treeHeight = df_dict[seg].treeHeight
                # ===== Draw connecting lines =====
                # 1st straight line, from tip.x to tre.TreeHt
                ax.plot([x0, x_offset+treeHeight], [y0, y0], 
                        c=cdict[clade_label], 
                        alpha=dotted_line_alpha,
                        linewidth=dotted_line_width)
                
                # can find tipname in the next tree
                if x1 != -1 and y1 != -1:
                    x1=x1+next_x_offset
                    y1=y1+next_y_offset
                    # Bent line, from tree.TreeHt to next_tree.treeHt
                    ax.plot([x_offset+treeHeight, next_x_offset], [y0, y1/y_zoom_dict[next_seg_dict[seg]]], 
                            c=cdict[clade_label], 
                            alpha=dotted_line_alpha,
                            linewidth=dotted_line_width)
                    # 2nd straight line
                    ax.plot([next_x_offset, x1], [y1/y_zoom_dict[next_seg_dict[seg]], y1/y_zoom_dict[next_seg_dict[seg]]], 
                            c=cdict[clade_label], 
                            alpha=dotted_line_alpha,
                            linewidth=dotted_line_width)
                

        elif isinstance(k,bt.node) or k.branchType=='node':
            line = np.array([[x+x_offset, k.children[0].y/y_zoom_dict[seg]+y_offset], [x+x_offset, k.children[-1].y/y_zoom_dict[seg]+y_offset]])
            lines_ls.append(line)

        line = np.array([[xp+x_offset, y+y_offset], [x+x_offset, y+y_offset]])
        lines_ls.append(line)


patch_collection = PatchCollection(patches_ls, color="k", zorder=10)
line_collection = LineCollection(lines_ls, lw=branchWidth,color='k', zorder=10)
ax.add_collection(patch_collection)
ax.add_collection(line_collection)


# ==================== Legend ====================

x0 = 0
x1 = 0.005
x_text = x0+0.007
y0 = 825
y_space = 65
# Specify clades in the legend, drawn from cdict.keys()
# To control order of iteration
clades_ls = ['3a1','3a2','3a2 (global)','3a','V1A.3', 'V1A', 'unassigned']


for i in range(len(clades_ls)):
    y = y0 - i*y_space
    ax.plot([x0, x1], [y, y], c=cdict[clades_ls[i]], lw=branchWidth*6)
    ax.text(x_text, y, clades_ls[i], 
            color="k", 
            verticalalignment="center", 
            fontsize=15)

# scale bar
y_bar = 825
ax.plot([x0+0.021, x1+0.021], [y_bar, y_bar], c="k", lw=branchWidth*6)
ax.text(x_text+0.021, y_bar, "0.005", 
        color="k", 
        verticalalignment="center", 
        fontsize=15)

ax.set_ylim(0, global_ySpan)
ax.set_xlim(0,cumulative_x_offset-0.005)
ax.set_xticks([])
ax.set_yticks([])
plt.axis("off")

print("Done in %.2fs" % (time.time() - t0))

plt.tight_layout()
plt.savefig("../results/tanglegram.pdf", bbox="inches")

plt.show()