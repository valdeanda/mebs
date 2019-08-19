
# coding: utf-8

#Import the libraries
import sys 
import matplotlib
import numpy as np
import matplotlib.pylab as plt
import pandas as pd
import seaborn as sns
from sys import argv

args = sys.argv
filename = args[1]
df = pd.read_table(filename)

#df=pd.read_table("../data/mats2plot.tab")
#df=pd.read_table("/home/val/src/metagenome_Pfam_score/scores2plot.tab")
df=df.rename(columns={ df.columns[0]: "met" })
met = df['met']
df.drop(labels=['met'], axis=1,inplace = True)
df.insert(5, 'met', met)
df=df.sort_values('sulfur', ascending=False)

df.columns[:-1]

# Make the PairGrid
g = sns.PairGrid(df.sort_values(["sulfur","carbon","nitrogen","iron","oxygen"], ascending=[False,False,False, False,False]),
                x_vars=df.columns[:-1], y_vars=['met'],
                size=10, aspect=.25)

# Draw a dot plot using the stripplot function
g.map(sns.stripplot, size=10, orient="h",
       palette="GnBu_d",edgecolor="gray")

# Use the same x axis limits on all columns and add better labels

g.set(xlabel="Score", ylabel="")

# Use semantically meaningful titles for the columns
titles = ["Sulfur", "Carbon", "Oxygen", "Iron", "Nitrogen"]

for ax, title in zip(g.axes.flat, titles):

    # Set a different title for each axes
    ax.set(title=title)
    

    # Make the grid horizontal instead of vertical
    ax.xaxis.grid(False)
    ax.yaxis.grid(True)

sns.despine(left=True, bottom=True)
plt.tight_layout()
plt.savefig(argv[1]+".png")
plt.close()






