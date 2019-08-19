# -*- coding: utf-8 -*-
# requires matplotlib >= v1.4

from sys import argv
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

if __name__ == '__main__':
    fname = argv[1]
    # Use . as decimal indicator
    data = pd.read_table(fname, index_col=0, na_values=['NA', "SIN DATO"], decimal='.')
    dataT = data.T

    # compute percentiles across Pfam families
    perc5arr = []
    perc95arr = []
    for col in range(0,dataT.columns.size-1):
      perc5 = np.nanpercentile(data.values[col],5)
      perc95 = np.nanpercentile(data.values[col],95)
      perc5arr.append(perc5)
      perc95arr.append(perc95)
      print (dataT.columns[col], ' 5-percentile=', perc5,' 95-percentile=', perc95)
   
    minperc5 = np.nanmin(perc5arr)
    maxperc95= np.nanmax(perc95arr) 
    print ('# min  5-percentile=', minperc5)
    print ('# max 95-percentile=', maxperc95)

    # plot boxplot of profiles with 5- and 95-percentiles whiskers
    plt.figure(figsize=(7,15))
    dataT.boxplot(grid=False, vert=False, whis=[5,95] )
    plt.axvline(0, alpha=0.5)
    plt.axvline(minperc5, alpha=0.5, color='black')
    plt.axvline(maxperc95, alpha=0.5, color='black')    
    plt.title("black lines: min 5% and max 95% percentiles")
    plt.yticks(size='xx-small')
    plt.xlabel('Entropy (bits)')
    plt.savefig(argv[1]+"_prof_box.png")
    plt.close()
 
    # Plot scatter plot
    means = data.mean(1)
    stds = data.std(1)
    df = pd.DataFrame([means, stds], index=['mean', 'std'])
    df = df.T
    df.plot(x='mean', y = 'std', kind='scatter')
    plt.tight_layout()
    plt.savefig(argv[1]+"_scatter.png")
    plt.close()
    
    
    


