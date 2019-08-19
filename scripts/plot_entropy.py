 # -*- coding: utf-8 -*-

from sys import argv
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
#sudo apt-get install python3-cairocffi
def round5(n):
    """Up round a decimal number to next, 0.5 multiple
    
    Arguments:
    - `n`: float
    """
    res = n % 0.5
    plus = 0.5 - res
    if n < 0:
        rounded = n - res
    elif n > 0:
        rounded = n + plus
    else:
        rounded = n
    return rounded

def plot_graph(dataframe):
    """Plot the data in  horizontal bars. Each column.
    
    Arguments:
    - `dataframe`: pandas data frame
    """
    # dimmentions in cm
    w = 12/2.5*2
    h = 15/2.5*2
    fig = plt.figure(figsize=(w, h))
    #
    _max = dataframe.max().max()
    _min = dataframe.min().min()
    elements=len(dataframe)
    columns = len(dataframe.columns)
    colors = plt.cm.rainbow(np.linspace(0, 1, columns))
    ploti = 101 + columns*10
    width = 1
    bottom = np.arange(0, elements, 1)
    for column in dataframe:
        data = dataframe[column]
        ax = plt.subplot(ploti)
        c = colors[(ploti % 10) - 1]
        ax.barh(bottom, data, height=width, lw=0.2, color=c, label=column)
        ax.set_ylim((0, elements))
        xmax = round5(_max)
        xmin = round5(_min)
        _xticks = np.arange(xmin, xmax+0.1, 0.5)
        plt.xticks(_xticks-0.05, _xticks, rotation=90, size=10)
        if ploti % 10 == 1:
            plt.yticks(np.arange(elements) + 0.5,
                       dataframe.index, size=6)
            plt.ylabel('Domain')
        elif ploti % 10 == columns/2:
            plt.xlabel('Entropy')
            ax.set_yticklabels([])
        else:
            ax.set_yticklabels([])
        # - title
        plt.title(column)
        # - subplot update
        ploti += 1
    fig.tight_layout()
    fig.subplots_adjust(wspace=0.18)        
    
def hmap(dataframe):
    """
    
    Arguments:
    - `datafreame`:
    """
    w = 7/2.5*2
    h = 15/2.5*2
    fig = plt.figure(figsize=(w, h))
    #
    _max = dataframe.max().max()
    _min = dataframe.min().min()
    elements=len(dataframe)
    columns = len(dataframe.columns)
    colors = plt.cm.rainbow(np.linspace(0, 1, columns))
    ploti = 101 + columns*10
    width = 1
    bottom = np.arange(0, elements, 1)    
    elements=len(dataframe)
    columns = len(dataframe.columns)
    _max = dataframe.max().max()
    cmap = plt.cm.seismic
    cmap.set_bad('g')
    # mask for nan values
    masked = np.ma.array(dataframe, mask=np.isnan(dataframe))
    plt.pcolor(masked, cmap=cmap, vmin=-_max,
                    vmax=_max, edgecolor='k')
    ax = plt.gca()
    plt.ylim((0, elements))
    plt.xlim((0, columns))
    plt.colorbar(fraction=0.05, shrink=.30)
    plt.xticks(np.arange(columns)+0.5, dataframe.columns,
               rotation=90)
    ax.xaxis.set_ticks_position('top')
    plt.yticks(np.arange(elements)+0.5, dataframe.index,
               size=6)
    plt.xlabel("Entropy")
    plt.ylabel("Domains")
    fig.subplots_adjust(bottom=0.05, top=0.95)
    
if __name__ == '__main__':
    if len(argv) >3:
      fname = argv[1]
      perc5rnd = float(argv[2])
      perc95rnd = float(argv[3])
    elif len(argv) > 1: 
      fname = argv[1]
      perc5rnd = 0.0
      perc95rnd = 0.0
    else:
      print ("""Usage:
        $python3 plot_entropy.py data.tab [random perc5] [random perc95]
        """)
      exit()

    # Use . as decimal indicator
    data = pd.read_table(fname, index_col=0, na_values=['NA', "SIN DATO"], decimal='.')
    # sort data by 'real' column
    #old version pandas  0.17 
    #data.sort('real', inplace=True)
    data.sort_values('real',inplace=True)
    plot_graph(data)
    plt.savefig(argv[1]+"_bar.png")
    plt.close()
    hmap(data)
    plt.savefig(argv[1]+"_hmap.png")
    plt.savefig(argv[1]+"_hmap.pdf")
    plt.close()
    # plot entropies histograms
    data.hist(figsize=(10, 10)) #,range=[-0.3, 1.3])
    plt.savefig(argv[1]+"_entropy_hist.png")
    plt.close()
    # plot boxplot of profiles
    #data.T.boxplot(figsize=(10,15), grid=False, vert=False )
    plt.figure(figsize=(7,15))

    # default boxplot with whiskers around 1.5x IQR
    data.T.boxplot(grid=False, vert=False )
    plt.axvline(0, alpha=0.5)

    if(perc95rnd < 0.0 or perc95rnd > 0.0):
      plt.axvline(perc5rnd, alpha=0.5, color='black')
      plt.axvline(perc95rnd, alpha=0.5, color='black') 

    plt.yticks(size='xx-small')
    plt.xlabel('Entropy (bits)')
    #plt.tight_layout()
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

    # Plot differential plot
    dif = data.T - data['real']
    dif = dif.drop('real') 
    meandif = dif.mean(1)
    stddif = dif.std(1)

    plt.errorbar(meandif.index, meandif, yerr=stddif, fmt='k--o')
    plt.xlabel('MSL', weight='bold')
    plt.ylabel('entropy difference with respect to real (bits)', weight='bold')
    plt.tight_layout()
    plt.savefig(argv[1]+"_differential.png")
    plt.close()















   
    
    


