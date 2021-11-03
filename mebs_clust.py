#!/usr/bin/env python3
 # -*- coding: utf-8 -*-
 #
 # ------------------------------
 # Name:     mebs_clust.py
 # Purpose:   Compute clustering analysis based on  presence abscence  matrix
 #
 # Authors:     acph - dragopoot@gmail.com and
 #              vydat - valdeanda@ciencias.unam.mx
 # Created:     2020
 # Licence:     GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
 # ----

# Import libraries

import argparse
from pathlib import Path
from argparse import RawDescriptionHelpFormatter
import numpy as np
import matplotlib.pylab as plt
from matplotlib import cm
import matplotlib as mpl
from matplotlib.lines import Line2D
import pandas as pd
from scipy import spatial
# from scipy.spatial.distance import pdist
# from scipy import cluster
from scipy.cluster import hierarchy

# Arguments and options
epilog = """Example:
$  python3 mebs_clust.py  pfamprofile.tsv"""

parser = argparse.ArgumentParser(description=__doc__, epilog=epilog,
                                 formatter_class=RawDescriptionHelpFormatter)

parser.add_argument('filename',
                    help="Input file pfam profile.")

parser.add_argument('--outdir', '-o', type=str, default=None,
                    help=('''Output folder [<filename>_mebs_clust]'''))

parser.add_argument('-im_format', '-f', default='pdf', type=str,
                    choices=['png', 'pdf', 'ps', 'eps', 'svg', 'tif',
                             'jpg'],
                    help='''Output format for images [pdf].''')

parser.add_argument('--im_res', '-r', default=300, type=int,
                    help='''Output resolution for images in
                    dot per inch (dpi) [dpi].''',
                    metavar='dpi')

parser.add_argument('--cutoff', '-co', default=0.5, type=float,
                    help=('''Cutoff threashold default 0.5'''))

parser.add_argument('--method', '-m',
                    choices=['single', 'complete', 'average', 'weighted',
                             'centroid', 'median', 'ward'],
                    default='ward',
                    type=str,
                    help="""methods for calculating the distance between
                    the newly formed clusters
                    [ward].""")

parser.add_argument('--distance', '-d',
                    choices=['braycurtis', 'canberra', 'chebyshev',
                             'cityblock', 'correlation',
                             'euclidean', 'jaccard', 'mahalanobis'],
                    default='jaccard',
                    type=str,
                    help="""The distance metric to use default
                    [jaccard].""")

parser.add_argument('--nolegend', '-nl', action='store_false',
                    default=True,
                    help="""If specified, the legend is not shown""")


args = parser.parse_args()

# END options

# Output dir name
if args.outdir == None:
    outpath_name = args.filename + '_mebs_clust'
    print('[OUTPUT]   Output not specified.')
else:
    outpath_name = args.outdir
print('[OUTPUT]   Storing result in : {}'.format(outpath_name))
outpath = Path(outpath_name)
if not outpath.exists():
    print('[OUTPUT]   ...', str(outpath), '-> does not exists: creating')
    outpath.mkdir()

# Variables
filename = args.filename
method = args.method
distance = args.distance
legend = args.nolegend

# Input files
df = pd.read_table(filename, index_col=0)
# convert to boolean
kobool = df > 0
kojacc = spatial.distance.pdist(kobool, metric=distance)
Z = hierarchy.linkage(kojacc, method=method, optimal_ordering=True)
maxdist = np.max(Z[:, 2])

# clustering (segmentación)
clust = hierarchy.fcluster(Z, maxdist*args.cutoff, criterion='distance')
n_clusts = len(np.unique(clust))       # number of clusters

# dendograms
# colors for clusters
# TODO: Add an option for selecting colormap
dend_colors = cm.jet(np.linspace(0, 1, n_clusts))
hexcolors = [mpl.colors.rgb2hex(rgb[:3]) for rgb in dend_colors]
hierarchy.set_link_color_palette(hexcolors)

# Plot
plt.figure(figsize=(15, 7))
dend = hierarchy.dendrogram(Z, color_threshold=maxdist*args.cutoff,
                            no_labels=True)
plt.axhline(maxdist*args.cutoff, ls='--', alpha=0.3, c='k')
plt.tight_layout()

# legend and colors
if legend:
    legend_elements = []
    for i, rgb in enumerate(dend_colors):
        label = f"Cluster {i+1}"
        element = Line2D([0], [0], color=rgb, label=label, lw=3)
        legend_elements.append(element)
    plt.legend(handles=legend_elements, loc=1)

print(f'[INFO] Output path: {outpath}')
outfname = outpath/(args.filename + f'.cutoff_{str(args.cutoff)}_' + "mebs_clust." + args.im_format)
plt.savefig(outfname,
            dpi=args.im_res, bbox_inches='tight')

# Store genomes of each of the  clusters
for c in np.unique(clust):
    print(f'Storing cluster {c} ...', end='\r')
    members = df[clust == c].index
    outfn = outpath/(args.filename + f'.cutoff_{str(args.cutoff)}_cluster{c}.txt')
    with open(outfn, 'w') as outf:
        for mem in members:
            outf.write(mem+'\n')
print("[END] mebs_cust.py done 9.......................\n"
      "[END] Please check your output files in folder '{}' :\n".format(str(outpath)))
