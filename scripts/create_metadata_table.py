#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
#------------------------------
# Name:
# Purpose:
#
# @uthor:   acph - dragopoot@gmail.com
#
# Created:     Wed Sep 23 13:42:39 CDT 2015
# Copyright:   (c) acph 2015
# Licence:     GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
#------------------------------
"""Creates a tabular and a pickle file that contains a table of metagenome 
metadata using the json mg-rast files contained in the specified directory."""

import argparse
import os
import json
import pandas as pd


def aparser():
    """Arguments parser function. Returns argumens object

    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directorypath',
                        help="""Directory that contains only the mg-rast json files. 
Each file corresponds to one metagenome.""")
    parser.add_argument('--excel', action='store_true',
                        help="""Creates an excel file of the data.""")
    args = parser.parse_args()
    return args


def env_package_data(infol):
    """Create a set of unique parameters measured among the
    metagenomes stored in the infol.

    The metadata data of each metagenome is stored y json format.

    The data is from MGRAST

    Arguments:
    - `infol`: Path where the metagenomes metadata is stored
    """
    files = os.listdir(infol)
    attr = set([])              # list of attributes
    # - for control observation
    # nattr = []
    for file_ in files:
        fname = os.path.join(infol, file_)
        with open(fname) as inf:
            txt = inf.read()
        data = json.loads(txt)
        if 'env_package' not in data['metadata']:
            continue
        dattr = data['metadata']['env_package']['data'].keys()
        # - for control observation
        # nattr.append(len(dattr))
        for value in dattr:
            attr.add(value)
    return attr


def main():
    """
    """
    # Create two dataframes, one for the general data and other
    # for the attributes of env_package
    args = aparser()
    infol = args.directorypath
    files = os.listdir(infol)
    # -- General dataframe
    nmet = len(files)
    index = pd.Index(files, name='metagenome')
    columns = ['biome', 'feature', 'material', 'env_package',
               'location', 'latitude', 'longitude', 'depth']
    generaldf = pd.DataFrame(index=index, columns=columns)
    # -- env_package attributes dataframe
    attr = list(env_package_data(infol))
    attrdf = pd.DataFrame(index=index, columns=attr)
    print('Reading files in:', infol)
    # iteration
    for meta in files:
        fname = os.path.join(infol, meta)
        with open(fname) as inf:
            data = json.loads(inf.read())
        for col in columns:
            gendic = data['metadata']['sample']['data']
            if col in gendic:
                generaldf[col][meta] = gendic[col]
                if col is 'env_package':
                    attdic = data['metadata']['env_package']['data']
                    for att, val in attdic.items():
                        attrdf[att][meta] = val
    # end iteration
    print("Saving files...")
    # save pandas
    prefix = infol.split('/')[-1]
    generaldf.to_pickle(prefix + '_generaldata.pk')
    generaldf.to_csv(prefix + '_generaldata.tab', sep='\t')
    attrdf.to_pickle(prefix + '_attributes.pk')
    attrdf.to_csv(prefix + '_attributes.tab', sep='\t')
    # save to excel
    if args.excel:
        excel_file = pd.ExcelWriter(prefix + '_data.xls')
        generaldf.to_excel(excel_file, 'general data')
        nullgeneral = pd.isnull(generaldf).sum()
        nullgeneral = pd.DataFrame(nullgeneral,
                                   columns=['Count of missing data'])
        nullgeneral.to_excel(excel_file, 'general-null')
        # -
        attrdf.to_excel(excel_file, 'env_package attributes')
        nullatt = pd.isnull(attrdf).sum()
        nullatt = pd.DataFrame(nullatt,
                               columns=['Count of missing data'])
        nullatt.to_excel(excel_file, 'attributes-null')
        # -
        excel_file.save()


if __name__ == '__main__':
    main()
    print("Done!!! Have a nice day.")
