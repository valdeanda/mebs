# coding: utf-8

#Import the libraries
import sys
from sys import argv
args = sys.argv

import numpy as np
import seaborn as sns
import matplotlib.pylab as plt
import pandas as pd
import numpy as np
import matplotlib
#Customize background and pallete colors 
sns.set_style("darkgrid")





filename = args[1]
df = pd.read_table(filename, index_col=0, sep="\t")

df_comp=df.drop(['sulfur', 'carbon','oxygen','iron','nitrogen','<sulfur comp>', '<carbon comp>'],axis=1)

df_comp.rename(columns={'sulfur_1': 'Sulfite oxidation',
'sulfur_2':'Thiosulfate oxidation', 
'sulfur_3':'Tetrathionate oxidation', 
'sulfur_4':'Tetrathionate reduction', 
'sulfur_5':'Sulfate reduction DS',
'sulfur_6':'Elemental sulfur reduction', 
'sulfur_7':'Thiosulfate disproportion', 
'sulfur_8':'Carbon disulfide oxidation', 
'sulfur_9':'Alkanesulfonate degradation', 
'sulfur_10':'Sulfate reduction A',
'sulfur_11':'Sulfide oxidation', 
'sulfur_12':'Cysteate oxidation', 
'sulfur_13':'Dimethylsulfone oxidation', 
'sulfur_14':'Sulfoacetate oxidation',
'sulfur_15':'Sulfolactate oxidation', 
'sulfur_16':'DMS oxidation',
'sulfur_17':'DMSP oxidation',
'sulfur_18':'MTP oxidation', 
'sulfur_19':'Suloacetaldehyde oxidation',
'sulfur_20':'Elemental sulfur oxidation',
'sulfur_21':'Elemental sulfur disproportion', 
'sulfur_22':'Methanesulfonate oxidation', 
'sulfur_23':'Taurine oxidation', 
'sulfur_24':'DMS methanogenesis', 
'sulfur_25':'MTP methanogesis', 
'sulfur_26':'Methanethiol methanogenesis', 
'sulfur_27':'Homotaurine degradation', 
'sulfur_28':'SQDG biosynthesis',
'sulfur_29':'Marker genes', 
'carbon_1':'coenzyme B/coenzyme M regeneration I (methanophenazine-dependent)',
'carbon_2': 'Methane oxidation, methanotroph, methane => formaldehyde',
'carbon_3': 'methanogenesis energy conservation', 
'carbon_4': 'Methanogenesis, acetate => methane (M00357)', 
'carbon_5':'Methanogenesis, methylamine/dimethylamine/trimethylamine => methane',
'carbon_6':'Methanogenesis from dimethylsulfide/methanethiol/methlthiolpropanoate => methane',
'carbon_7':'Methanogenesis, CO2 => methane',
'carbon_8':'methanogenesis from acetate reductive acetyl coenzyme A pathway II (autotrophic methanogens)',
'carbon_9':'Methanogenesis, methanol => methane',
'carbon_10':'methylamine degradation',
'carbon_12':'methyl-coenzyme M oxidation to CO2',
'carbon_13':'methyl-coenzyme M reduction to methane' },inplace=True) 

df_comp=df_comp.T
sns.set(font_scale=1)
axs = sns.clustermap(df_comp, col_cluster=True, linewidths=0.6,cmap=sns.color_palette("RdBu_r", 100),
                     figsize=(15,12))
plt.tight_layout()
plt.savefig(argv[1]+".png", bbox_inches='tight', dpi=500)
plt.close()
