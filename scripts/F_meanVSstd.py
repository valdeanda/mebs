#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# ------------------------------
# Name:     F_meanVSstd.py
# Purpose:  Mean vs std figure of profiles Using ward linkage hyerarchical clustering
#
# @uthor:      acph - dragopoot@gmail.com
#
# Created:     2015
# Copyright:   (c) acph 2015
# Licence:     GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
# ------------------------------
"""Mean vs standard deviation figure of profiles and clustering. Creates a file for
each cluster that contains the list of profiles that are included."""

import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from sklearn import cluster   # , datasets
from sklearn.neighbors import kneighbors_graph
from sklearn.preprocessing import StandardScaler


def limits(array, percentage=0.01):
    """Computes plot limits to plot the data of an array.

    Parameters
    ----------
    array : 1D array
    percentage : Fraction of the array to add to the limits


    Returns
    -------
    out : Returns a 2 values tuple. First value is low limit and second
          value is the high limit

    """
    max_ = array.max()
    min_ = array.min()
    ad = (max_ - min_) * percentage
    if max_ == 0.0:
        max_ = 0 + ad
    if min_ == 0.0:
        min_ = 0 - ad
    low = min_ - ad
    high = max_ + ad
    return low, high

# options
epilog = """Example:

$ python3 F_meanVSstd entropies_matrix_entropies.tab -o figure.png"""

parser = argparse.ArgumentParser(description=__doc__, epilog=epilog)
parser.add_argument('filename',
                    help="Input file in tabular format. Rows are pfam families and " +
                    "columns are metagenome fragment (reads) length.")
parser.add_argument(
    '-o', '--out_fig', help='Stores the figure in the specified file (and format).')
parser.add_argument('--dpi', type=int, default=300,
                    help='Resolution for output figure file [300].')
parser.add_argument('-v', '--variation', default='std',
                    choices=['std', 'cv', 'id', 'range'],
                    help='Select the measurement of variation to plot in y axis [std]: standard devitation (std), coefficient of variation (cv), index of dispersion (id) or range. cv and id cannot be used in variables with negative values.')
parser.add_argument('-k', type=int, choices=range(2, 9), default=3,
                    help='Number of k-means clusters [3].')
parser.add_argument('--plot-random', default=None, metavar='DIRECTORY',
                    help='Folder where the *.tab files containing random samples are stored.')
parser.add_argument('-c', '--cluster-alg', choices=['ward', 'birch'],
                    default='ward',
                    help='Chose clustering algorithm [ward]. Ward linked hierarchical clustering or birch clustering.')
parser.add_argument('--labels', type=int, metavar='CLUSTER',
                    help="Plot the labels of the points in the specified cluster.")
args = parser.parse_args()

# input file
# fname = 'matrices_pfam_entropies.tab'
# fname = 'matrices_curadas_sep_entropies.tab'
fname = args.filename

data = pd.read_table(fname, index_col=0, na_values=['NA', "SIN DATO"])
#                     decimal='.')
means = np.array(data.mean(1))
stds = np.array(data.std(1))
mask1 = ~np.isnan(means)
mask2 = ~np.isnan(stds)
mask = mask1 * mask2
data = data.ix[mask]
means = np.array(data.mean(1))
stds = np.array(data.std(1))

#############
# Variation #
#############
# cv = stds / means                    # Coefficient of variation
# ID = stds**2 / means                 # Index of dispersion
#_range = np.array(data).ptp(1)       # Data Range
if args.variation == 'std':
    var_ = stds
    y_label = 'Entropy standard deviation'
elif args.variation == 'cv':
    var_ = stds / means
    y_label = 'Entropy coefficient of variation'
elif args.variation == 'id':
    var_ = stds**2 / means
    y_label = 'Entropy index of dispersion'
elif args.variation == 'range':
    var_ = np.array(data).ptp(1)
    y_label = 'Entropy range'

##############
# Clustering #
##############
x = np.vstack((means, var_)).T
k = args.k
# noramlize data
X = StandardScaler().fit_transform(x)

# connectivity for ward clustering: search for neighbors for each point
connectivity = kneighbors_graph(X, n_neighbors=10, include_self=False)

# ward clustering
ward = cluster.AgglomerativeClustering(
    n_clusters=k, linkage='ward', connectivity=connectivity)
# birch clustering
birch = cluster.Birch(n_clusters=k)

# algorithm selection
algorithms = {'ward': ward, 'birch': birch}
alg = algorithms[args.cluster_alg]
# fit
alg.fit(X)
# cluster labels
y_pred = alg.labels_
clusts = np.unique(y_pred)
cs_ = ['b', 'g', 'r', 'y', 'k', 'c', 'm', 'grey']

# figure
fig = plt.figure()
# Ax positions
# scat_pos = [0.15, 0.15, 0.7, 0.7]
# xbox_pos = [0.15, 0.85, 0.7, 0.1]
# ybox_pos = [0.85, 0.15, 0.1, 0.7]


axscatter = fig.add_subplot(111)
# axscatter = fig.add_axes(scat_pos, frameon=True)
# axxbox = fig.add_axes(xbox_pos, frameon=False)
# axybox = fig.add_axes(ybox_pos, frameon=False)

# plot random
if args.plot_random:
    dataframes = {}
    infolder = args.plot_random
    files = os.listdir(infolder)
    files = [f for f in files if '.tab' in f]
    for f in files:
        path = os.path.join(infolder, f)
        key = f.split('.')[-2]
        df = pd.read_table(path, index_col=0, na_values=['NA'])
        del df[df.columns[0]]
        dataframes[key] = df
    panel = pd.Panel(dataframes)
    r_means = panel.mean(0)
    r_stds = panel.std(0)
    if args.variation == 'std':
        r_var_ = r_stds
    elif args.variation == 'cv':
        r_var_ = r_stds / r_means
    elif args.variation == 'id':
        r_var_ = r_stds**2 / r_means
    elif args.variation == 'range':
        r_var_ = panel.max(0) - panel.min(0)
#    import seaborn as sns
#    sns.kdeplot(r_means, r_variation, ax=axscatter)
    if args.variation in ['cv', 'id']:
        axscatter.scatter(r_means, r_var_, color='mistyrose',
                          label='Random samples', alpha=1)
        print("Warning!!! - The 2d histogram of random data only",
              "works with standard deviation and range.\n",
              "Plotting simple scatter plot for random data.")

    else:
        from matplotlib.colors import LogNorm
        h2d = axscatter.hist2d(r_means.get_values().flatten(),
                               r_var_.get_values().flatten(),
                               bins=50, norm=LogNorm(), cmap='RdPu', alpha=0.7,
                               label='Random samples')
        mappable = h2d[3]
        cbar = plt.colorbar(mappable, pad=0.02, fraction=0.1, use_gridspec=False,
                            anchor=(0.0, 0.0),
                            shrink=0.5)
        cbar.set_label('Random sample frequency (log)', fontsize='small')

# scatter plot
for i in clusts:
    mask = y_pred == i
    data_ = x[mask]
    axscatter.scatter(data_[:, 0], data_[:, 1], alpha=0.5,
                      color=cs_[i], label='Cluster {}'.format(i))

leg = axscatter.legend(fontsize='small')
leg.draggable()
ylims = limits(x[:, 1])
xlims = limits(x[:, 0])
axscatter.set_ylim(ylims)
axscatter.set_xlim(xlims)

if type(args.labels) == int:
    clust2 = x[y_pred == args.labels]
    labels = data.index[y_pred == args.labels]
    for i in range(len(labels)):
        axscatter.annotate(labels[i], clust2[i],
                           alpha=0.5, fontsize='xx-small')


axscatter.set_xlabel("Entropy mean", fontweight='bold')
axscatter.set_ylabel(y_label, fontweight='bold')

# box plots
# bpx = axxbox.boxplot(x[:, 0], vert=False)
# axxbox.set_xticks([])
# axxbox.set_yticks([])
# axxbox.set_xlim(xlims)
# axxbox.set_ylim(0.9, 1.1)

# bpy = axybox.boxplot(x[:, 1], vert=True)
# axybox.set_xticks([])
# axybox.set_yticks([])
# axybox.set_ylim(ylims)
# axybox.set_xlim(0.9, 1.1)

# plt.setp(bpx['boxes'], color='black', linewidth=1.5)
# plt.setp(bpx['whiskers'], color='black', linewidth=1.5)
# plt.setp(bpx['caps'], color='black', linewidth=1.5)
# plt.setp(bpx['fliers'], color='black')

# plt.setp(bpy['boxes'], color='black', linewidth=1.5)
# plt.setp(bpy['whiskers'], color='black', linewidth=1.5)
# plt.setp(bpy['caps'], color='black', linewidth=1.5)
# plt.setp(bpy['fliers'], color='black')


# save clusters
for i in clusts:
    mask = y_pred == i
    dset = data.ix[mask]
    m = dset.mean(1)
    s = dset.std(1)
    _cv = s / m
    _id = s**2 / m
    _r = dset.max(1) - dset.min(1)  # range
    df = pd.concat((m, s, _cv, _id, _r), 1,
                   keys=['mean', 'std', 'cv', 'id', 'range'])
    fname = 'cluster_{}_{}_k{}_{}.tab'.format(i, args.variation,
                                              args.k, args.cluster_alg)
    df.to_csv(fname, sep='\t')

if args.out_fig:
    plt.savefig(args.out_fig, dpi=args.dpi)
else:
    plt.show()
