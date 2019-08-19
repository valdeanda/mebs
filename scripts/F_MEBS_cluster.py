#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# ------------------------------
# Name:     F_MEBS_cluster.py
# Purpose:Creates a figure with the low dimention projection and clustering of the data
#         in a MEBS result file
#
# @uthor:      acootph - dragopoot@gmail.com
#
# Created:     2018
# Copyright:   (c) acph 2018
# Licence:     GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
# ------------------------------

"""Creates a figure with the low dimention projection and clustering of the data
in a MEBS result file. The clustering is computed from the stadardized and scaled
original data, not with the dimensionaly reduced."""

# import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from sklearn import cluster
from sklearn.preprocessing import StandardScaler, MaxAbsScaler, RobustScaler
from sklearn.manifold import TSNE, Isomap
from sklearn.decomposition import PCA

global projection
global clusterig

projection = {'tsne': TSNE(),
              'isomap': Isomap(),
              'pca': PCA(),
              }

clustering = {'ward': cluster.AgglomerativeClustering(),
              'kmeans': cluster.KMeans(),
              'spectral': cluster.SpectralClustering(),
              'meanshift': cluster.MeanShift(),
              }

scalers = {'std': StandardScaler(),
           'robust': RobustScaler(),
           'max': MaxAbsScaler()}


def proj_clust(df, k, n_components=2, scale='std', proj='tsne', clust='ward'):
    """Projecton and clustering of the given DataFrame.

    Parameters
    ----------
    df : pandas.DataFrame
        Pandas numerical DataFrame.
    k : int or float
        Cluster parameter. Must be int if clust is 'ward', 'kmeans' or
        'spectral' and correspond to the desired number of clusters. Must
        be float if the clust is 'meanshift' and represent the quantile
        value for bandwidth estimation. Ultimately, affects the number
        of clusters identified.
    n_components : int, [2|3]
        Number of dimensions to reduce the data.
    scale: str
        Name of the sacling method. Values specified in scalers
        dictionary.
    proj : str
        Name of the projection method. Values specified in projection
        dictionary.
    clust : str
        Name of the clustering method. Values specified in clustering
        dictionary.

    Returns
    -------
    projection : ndarray
        An array that contain the low dimension projection of the original
        data.
    clust_model : Object, sklearn model
        Trained clustering model.
    """
    # scaling data
    if scale == 'none':
        X = df.values
    else:
        scaler = scalers[scale]
        X = scaler.fit_transform(df.values)
    # projection
    projm = projection[proj]
    projm.n_components = n_components
    X_proj = projm.fit_transform(X)
    # Cluster
    if clust != 'meanshift':
        clustm = clustering[clust]
        clustm.n_clusters = k
    else:
        assert type(k) == float, 'meanshift k parameter must be a float'
        bandwidth = cluster.estimate_bandwidth(X, quantile=k)
        clustm = clustering[clust]
        clustm.bandwidth = bandwidth
    clustm.fit(X)
    return X_proj, clustm


def plot_projclust(X, clust_model, filename,
                   xlabel='Dimension 1',
                   ylabel='Dimension 2',
                   zlabel='Dimension 3',
                   leg_title='Model',
                   dpi=300):
    """Plot the projection and clustering of the data.

    Parameters
    ----------
    X : ndarray
        ndarray with the low dimension projection.
    clust_model : Object, sklearn model
        Trained clustering model.
    filename : str
        String that contain the path to a file to save
        the plot. Must end in a image extention, like .png.
    dpi : int > 0
        The resolution in dots per inch of the plot.
    xlabel, ylabel ,zlabel: str
        Names for x, y and z axis respectively.
    leg_title : str
        Name for the legend title.
    """
    # preparing
    labels = clust_model.labels_
    clusts = np.unique(clust_model.labels_)
    fig = plt.figure()
    # creating image
    if X.shape[1] == 2:
        ax = fig.add_subplot(111)
        for c in clusts:
            mask = labels == c
            X_c = X[mask]
            label = 'Cluster {}'.format(c)
            ax.scatter(*X_c.T, s=6, alpha=0.3,
                       label=label)
    else:
        assert X.shape[1] in [2, 3], 'Projection dimensions must be 2 or 3'
        print('not yet implemented')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if X.shape[1] == 3:
        ax.set_zlabel(zlabel)
    ax.legend(title=leg_title, loc=1)
    # save figure
    plt.savefig(filename, dpi=dpi)


def plot_all(df, fname, dpi=300, scale='std', k=4, kb=0.2):
    """Make a plot with all the projections and clustering algorithms allowed
    :p"""
    # preparing
    n = len(projection)
    m = len(clustering)
    lproj = list(projection.keys())
    lclust = list(clustering.keys())
    fig = plt.figure()
    count = 1
    ks = {'ward': k,
          'kmeans': k,
          'spectral': k,
          'meanshift': kb}
    # creating
    for i in range(n):
        p_ = lproj[i]
        for j in range(m):
            print('[PlotAll] Working ... {}/{}'.format(count, n * m),
                  flush=True, end='\r')
            c_ = lclust[j]
            ax = fig.add_subplot(n, m, count)
            # default k accepting a list form command line?
            k_ = ks[c_]
            # processing
            X_p, cmodel = proj_clust(df, k=k_,
                                     scale=scale,
                                     proj=p_, clust=c_)
            labels = cmodel.labels_
            clab = np.unique(cmodel.labels_)
            # creating image
            for cl in clab:
                mask = labels == cl
                X_c = X_p[mask]
                label = 'Cluster {}'.format(cl)
                ax.scatter(*X_c.T, s=6, alpha=0.3,
                           label=label)
                plt.title(scale + '+' + p_ + ' + ' + c_,
                          fontsize='xx-small')
            # plt.legend()
            plt.xticks(fontsize='xx-small')
            plt.yticks(fontsize='xx-small')
            # end processing
            count += 1
    plt.tight_layout()
    plt.savefig(fname, dpi=dpi)


def save_clust(df, clustm, basename):
    """Save the clusters generated from a clustering analisis
    in separeate files.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame with the original data. Index must be the labels
        of the observations.
    clustm : sklearn model
        Sklearn trainded clustering model.
    basename : str
        Base name for the files with the clusters


    Returns
    -------
    out : None

    """
    clusters_ = np.unique(clustm.labels_)
    samples = df.index.values
    for i in clusters_:
        mask = clustm.labels_ == i
        subsamp = samples[mask]
        fname = basename + "_cluster{}.txt".format(i)
        with open(fname, 'w') as outf:
            for samp in subsamp:
                line = samp + '\n'
                outf.write(line)


def arguments():
    "Program arguments"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('filename',
                        help="""Input file, MEBS result derived file.
TAB separated where the irst column is the sample name and the first
row  is the columns names. Each column corresponds to the MEBS score
for each pathway. The data is standardized and scaled before processing""")
    parser.add_argument('-k', '--kparam', default='4',
                        help="""Parameter for controling the number of clusters.
                        For ward, kmeans and spectral algortithms must be a
                        small integer that determines the exact number of
                        clusters. For meanshift,it must be a floating point
                        number in the range 0-1 that determines de variability
                        allowed to compute the clusters. Default 0.2 for
                        meanshift and 4 for the other algorithms.""")
    parser.add_argument('-s', '--scaler',
                        choices=['std', 'robust', 'max', 'none'],
                        default='std',
                        help="""Method for scale/standardize the data.
                        Standardization is usually needed for machine
                        learning methods. Based on scikit-learn, preprocessing
                        methods.
                        NOTE: 'none' option uses the
                        raw data for the analysis [std].""")
    parser.add_argument('-p', '--projection',
                        choices=['tsne', 'pca', 'isomap'],
                        default='tsne',
                        help="""Method for the projection for the data in 2 dimension
                        [tsne].""")
    parser.add_argument('-c', '--clustering',
                        choices=['ward', 'kmeans',
                                 'meanshift', 'spectral'],
                        default='ward',
                        help="""Method for data clustering
                        [ward].""")
    parser.add_argument('--all', action='store_true',
                        help="""All analysis in one figure using the default values
                        k=4 for ward, kmeans and spectral clustering and k=0.2
                        for meanshift.
                        This option allows to visualize the geneal behavior of
                        the methods for the subsequent selection of the desired
                        value. Does not store the clusters. May take some time :P""")
    parser.add_argument('--seed', type=int, default=2332,
                        help='''When applicable,
                        random number seed for all the analysis.''')
    parser.add_argument('-n', '--n-components', type=int, choices=[2, 3],
                        default=2,
                        help='''Number of dimensions for dimensionality
                        reduction. WARNING: just now the only valid option
                        is 2''')
    parser.add_argument('--dpi', type=int, default=300,
                        help='''Resolution of images in dot per inch [300].''')
    parser.add_argument('--im_format', '-f', default='png', type=str,
                        choices=['png', 'pdf', 'ps', 'eps',
                                 'svg', 'tif', 'jpg'],
                        help='''Output format for images [png].''')

    args = parser.parse_args()
    return args


def config_models(args):
    """General model configurartion of the  models using command line options
    njobs, random_state, etc

    Parameters
    ----------
    args : argparse.Namespace

    Returns
    -------
    out : None,
        global variables configured

    """
    projection = {'tsne': TSNE(random_state=args.seed),
                  'isomap': Isomap(n_jobs=-1),
                  'pca': PCA(random_state=args.seed),
                  }

    clustering = {'ward': cluster.AgglomerativeClustering(),
                  'kmeans': cluster.KMeans(random_state=args.seed,
                                           n_jobs=-1),
                  'spectral': cluster.SpectralClustering(random_state=args.seed,
                                                         n_jobs=-1),
                  'meanshift': cluster.MeanShift(n_jobs=-1),
                  }


def main():
    args = arguments()
    # "temporal" fix to meanshift default value
    if args.clustering == 'meanshift':
        args.kparam = '0.2'
    # ends temporal fix (lol)
    if '.' in args.kparam:
        args.kparam = float(args.kparam)
    else:
        args.kparam = int(args.kparam)
    df = pd.read_table(args.filename, index_col=0)
    basename = args.filename.split('.')[0]
    config_models(args)
    if args.all:
        print('[PlotAll] Creating all models...')
        allname = basename + '_plot_all_{}.{}'.format(args.scaler,
                                                      args.im_format)
        plot_all(df, allname, dpi=args.dpi,
                 scale=args.scaler,
                 k=args.kparam, kb=0.2)
        print('[PlotAll] Plotting and saving in {}...'.format(allname))
    else:
        # creating models
        print('[RedClust] Creating models...')
        assert args.n_components == 2, '3D projections still not implemented'
        proj_, cm_ = proj_clust(df, args.kparam,
                                n_components=args.n_components,
                                scale=args.scaler,
                                proj=args.projection,
                                clust=args.clustering)
        # creatig image
        fname = basename + '_{}_{}_{}.{}'.format(args.scaler,
                                                 args.clustering,
                                                 args.projection,
                                                 args.im_format)
        print('[RedClust] Plotting in : {}'.format(fname))
        axislabeldic = {'tsne': 'tSNE ',
                        'pca': 'PCA ',
                        'isomap': 'Isomap'}
        lab = axislabeldic[args.projection]
        plot_projclust(proj_, cm_, fname, dpi=args.dpi,
                       xlabel=lab + '1',
                       ylabel=lab + '2',
                       zlabel=lab + '3',
                       leg_title=args.clustering)
        # saving clusters
        print('[RedClust] Saving clusters: ...')
        fname = basename + '_{}_{}'.format(args.clustering,
                                           args.projection)
        save_clust(df, cm_, fname)


if __name__ == '__main__':
    main()
    print('[DONE] Have a nice day')
