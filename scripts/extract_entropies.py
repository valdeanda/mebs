#!/usr/bin/python
# -*- coding: utf-8 -*-
#
#------------------------------
# Name:     extract_entropies.py
# Purpose:  Extract the entropies valies from domains profiles
#           files
# 
# @uthor:   acph, converted to python3 by val and bruno
#
# Created:     mié abr 22 15:01:44 CDT 2015
# Copyright:   (c) acph 2015
# Licence:     GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
#------------------------------
"""Extract the entropies valies from domains profiles
 files """

import os
import re
from sys import argv
import pandas as pd

def extract_entropies(fname):
    """Extract the profile list and the entropies list
    from a domains profile (presence/absence) file

    Return a list of families names and the corrspondint list
    of entropy values
    
    Arguments:
    - `fname`: Domain profile file
    """
    tabline = re.compile('\S*?(\t\S*?)+\n')
    with open(fname) as inf:
        flag = True
        for line in inf:
            if tabline.search(line) and flag == True:
                line = line.strip('\n')
                profiles = line.split('\t')
                flag = False
            elif 'rel_entropy' in line:
                line = line.strip('\n')
                entropies = line.split('\t')
                break
    profiles = profiles[1:]
    entropies = entropies[1:]
    return profiles, entropies


def main():
    """
    """
    print("**********")
    print("Warning!!! All the files need to have the same")
    print("number of domains(profiles) in the same column order.")
    print("This script assumes that these considerations are")
    print("true, so it cannot find errors in the input files format")
    print("********** \n")
    infolder = argv[1]
    os.chdir(infolder)
    fnames = os.listdir('.')
    files = {}
    fnre = re.compile('_size([0-9]+)_')
    for fname in fnames:
        if 'size' not in fname:
            files['real'] = fname
        else:
            number = fnre.findall(fname)[0]
            files[number] = fname
            
    columns = ['real', '30', '60', '100', '150', '200',
               '250', '300']
    d = []
    for col in columns:
        fname = files[col]
        print("Working with {}...".format(fname))
        profiles, entropies = extract_entropies(fname)
        d.append(entropies)
    df = pd.DataFrame(d, index=columns, columns=profiles)
    df = df.T
    os.chdir('..')
    outfname = infolder.strip('/')
    outfname += '_entropies.tab'
    df.to_csv(outfname, sep='\t')
    

if __name__ == '__main__':
    if len(argv) != 2:
        print("""Usage:
        $python3 extract_entropies.py files_folder

        the files_folder must contain only the profiles files""")
    main()
    print("Done!! have a nice day")
    
