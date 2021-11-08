#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# ------------------------------
# Name:     mebs_vis.py
# Purpose:  Parse and compute some graphs derived from mebs.pl output
#
# Authors:     acph - dragopoot@gmail.com and vydat - valdeanda@ciencias.unam.mx
# Created:     2018
# Licence:     GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
# Last updated: October 2021
# Updates: Oct 2021: Genomic completeness based on marker genes module and separated file for normalized mebs scores
# ------------------------------

""" Parse mebs.pl output and creates several files and figures:
-File to map mebs normalized values to itol         => itol_mebs.txt
-File with the metabolic completeness with names    => input+completenes.tab
-File to be the output of F_MEBS_cluster.py -s none => input+2_cluster_mebs.txt
-Heatmap with normalized mebs values                => inputmebs_heatmap.png
-Heatmap with metabolic completness of S and C      => input+comp_heatmap.png
-Barplot with normalized mebs values                => input+barplot.png
-Genomic completeness based on marker genes         => input+genomic_completeness.tab
-Normalized mebs scores                             => input+norm_mebs.tab
"""

# Import libraries
import argparse
from pathlib import Path
from argparse import RawDescriptionHelpFormatter
# import numpy as np
import matplotlib.pylab as plt
import pandas as pd
import seaborn as sns
from subprocess import run

################################
# Script arguments and options #
################################
epilog = """Example:
$  python3 mebs_vis.py gen_test.tsv """

parser = argparse.ArgumentParser(description=__doc__, epilog=epilog,
								 formatter_class=RawDescriptionHelpFormatter)
parser.add_argument('filename',
					help="Input file derived from mebs.pl using -comp option.")
parser.add_argument('-o', '--outdir', type=str, default=None,
					help=('''Output folder [<filename>_mebs_vis]'''))

parser.add_argument('-im_format', '-f', default='png', type=str,
					choices=['png', 'pdf', 'ps', 'eps', 'svg', 'tif', 'jpg'],
					help='''Output format for images [png].''')
parser.add_argument('--im_res', '-r', default=300, type=int,
					help='''Output resolution for images in
					dot per inch (dpi) [dpi].''',
					metavar='dpi')


args = parser.parse_args()
# END options #


# Output dir name
if args.outdir == None:
	outpath_name = args.filename + '_mebs_vis'
	print('[OUTPUT]   Output not specified.')
else:
	outpath_name = args.outdir
print('[OUTPUT]   Storing result in : {}'.format(outpath_name))
outpath = Path(outpath_name)
if not outpath.exists():
	print('[OUTPUT]   ...', str(outpath), '-> does not exists: creating')
	outpath.mkdir()


############################################
# Remove asterisks from original mebs file #
############################################
noast_fname = outpath/( args.filename + '.noa')
cmd = "sed s/\*//g {} > {}".format(args.filename,
								   noast_fname)
run(cmd, shell=True)

# END remove asterisks #

df = pd.read_csv(noast_fname, index_col=0, sep='\t')

# Values obtained  by summing positive and negative entropies of each cycle
sval = [16.018, -6.527000000000001]
cval = [85.33199999999998, -6.555]
nval = [22.079000000000004, -0.7040000000000001]
feval = [10.463999999999999, -1.188]
oval = [10.702999999999998, -2.317]


# # Normalize from 0-1

def sulfur_per(x):
	if x >= 0:
		pct = float(x / sval[0])
		return round(pct, 2)
	elif x < 0:
		neg = float(abs(x / sval[1]))
		return round(neg, 2)


def carbon_per(x):
	if x >= 0:
		pct = float(x / cval[0])
		return round(pct, 2)
	elif x < 0:
		neg = float(abs(x / cval[1]))
		return round(neg, 2)


def nitrogen_per(x):
	if x >= 0:
		pct = float(x / nval[0])
		return round(pct, 2)
	elif x < 0:
		neg = float(abs(x / nval[1]))
		return round(neg, 2)


def iron_per(x):
	if x >= 0:
		pct = float(x / feval[0])
		return round(pct, 2)
	elif x < 0:
		neg = float(abs(x / feval[1]))
		return round(neg, 2)


def oxygen_per(x):
	if x >= 0:
		pct = float(x / oval[0])
		return round(pct, 2)
	elif x < 0:
		neg = float(abs(x / feval[1]))
		return round(neg, 2)

# Create a new column with normalized values


df['S'] = df.sulfur.apply(sulfur_per)
df['C'] = df.carbon.apply(carbon_per)
df['O'] = df.oxygen.apply(oxygen_per)
df['Fe'] = df.iron.apply(iron_per)
df['N'] = df.nitrogen.apply(nitrogen_per)

df_new = df[['S', 'C', 'O', 'Fe', 'N']]

# Send new dataframe only with normalized values
df_new.to_csv(outpath/(args.filename + "_norm_mebs.tab"), sep="\t")

# Create itol file with normalized values

outfilename = outpath/(args.filename + '_itol_mebs.txt')
infile = 'data2vis/dataset_heatmap_template.txt'
outfile = open(outfilename, 'w')

# Modify FIELD_LABELS
with open(infile) as inf:
	for line in inf:
		if 'FIELD_LABELS ' in line:
			new_line = ['FIELD_LABELS'] + list(df_new.columns)
			new_line = ' '.join(new_line) + '\n'   # spaces
			outfile.write(new_line)
		else:
			outfile.write(line)
# new line
outfile.write('\n')

# new df

for ind_ in df_new.index:
	l = [str(i) for i in df_new.loc[ind_]]
	line = ind_ + ' ' + ' '.join(l) + '\n'
	outfile.write(line)
outfile.close()


# Create file to be input of   F_MEBS_cluster.py  using -s none option

outfilename = outpath/(args.filename + '_2_cluster_mebs.txt')
infile = 'data2vis/mebs.gen.nr.norm.tab'
outfile = open(outfilename, 'w')

# Modify FIELD_LABELS
with open(infile) as inf:
	for line in inf:
		if 'nitrogen ' in line:
			new_line = ['nitrogen'] + list(df_new.columns)
			new_line = ' '.join(new_line) + '\t'   # separated by tab
			outfile.write(new_line)
		else:
			outfile.write(line)
# new line
# outfile.write('\n')

# new df

for ind_ in df_new.index:
	l = [str(i) for i in df_new.loc[ind_]]
	line = ind_ + '\t' + '\t'.join(l) + '\n'
	outfile.write(line)
outfile.close()


# Create the pfam metabolic completeness file

df_comp = df.drop(['sulfur', 'carbon', 'oxygen', 'iron', 'nitrogen',
				   '<sulfur comp>', '<carbon comp>', '<nitrogen comp>',
				   '<iron comp>', 'markers', '<markers comp>',
				   'S', 'C', 'N', 'O', 'Fe'], axis=1)
df_comp.rename(columns={'sulfur_1': 'aprAB(Marker_gene)',
						'sulfur_2': 'Apt/Sat(Marker_gene)',
						'sulfur_3': 'DsrABC(Marker_gene)',
						'sulfur_4': 'Sulfur_oxidation(Sox_system)',
						'sulfur_5': 'Sulfur_oxidation(Sor_system)',
						'sulfur_6': 'Sulfur_oxidation(FccB)',
						'sulfur_7': 'Sulfur_oxidation(DoxAD)',
						'sulfur_8': 'Sulfur_oxidation(DsrEFH)',
						'sulfur_9': 'DsrKMJOP(Marker_gene)',
						'sulfur_10': 'Sulfur_reduction(QmoABC)',
						'sulfur_11': 'Sulfur_oxidation(Puf_reaction_center)',
						'sulfur_12': 'Sulfur_assimilation(CysACDJNPQU)',
						'sulfur_13': 'Tetrathionate_reduction(asrABC)',
						'sulfur_14': 'Tetrathionate_reduction(ttrABC)',
						'sulfur_15': 'Thiosulfate_disproportionation(phsABC)',
						'sulfur_16': 'Thiosulfate_disproportionation(Rhodanase)',
						'sulfur_17': 'S°_reduction(hydACD)',
						'sulfur_18': 'S°_reduction(sreABC)',
						'sulfur_19': 'DMS_degradation(DdhABC)',
						'sulfur_20': 'DMS_degradation(DsoABCDEF)',
						'sulfur_21': 'DMS_degradation(DmoAB)',
						'sulfur_22': 'Sulfoacetaldehyde_degradation(isfD)',
						'sulfur_23': 'Sulfoacetaldehyde_degradation(Xsc)',
						'sulfur_24': 'Sulfoacetaldehyde_degradation(SafD)',
						'sulfur_25': 'Methanesulfonate_degradation',
						'sulfur_26': 'Sulfolactate_degradation',
						'sulfur_27': 'Taurine_degradation',
						'carbon_1': 'coB/coM_regeneration',
						'carbon_2': 'Methane_oxidation',
						'carbon_3': 'Methanogenesis',
						'carbon_4': 'Methanogenesis(methanol)',
						'carbon_5': 'Methylamine_degradation',
						'carbon_6': 'mcrABC(Marker_gene)',
						'nitrogen_1': 'Ammonia_assimilation_I',
						'nitrogen_2': 'Ammonia_assimilation_II',
						'nitrogen_3': 'L-glutamine_biosynthesis_I',
						'nitrogen_4': 'Superpathway_Ammonia_assimilation',
						'nitrogen_5': 'Ammonia_oxidation_I(aerobic)',
						'nitrogen_6': 'Ammonia_oxidation_II(anaerobic)',
						'nitrogen_7': 'Ammonia_oxidation_IV(Autotrophic_Ammonia_oxidizers)',
						'nitrogen_8': 'Nitrifier_denitrification',
						'nitrogen_9': 'Nitrate_reduction_I(denitrification)',
						'nitrogen_10': 'Nitrate_reduction_III(denitrification)',
						'nitrogen_11': 'Nitrate_reductionIV(dissimilatory)',
						'nitrogen_12': 'Nitrate_reductionV(assimilatory)',
						'nitrogen_13': 'Nitrate_reductionVI(assimilatory)',
						'nitrogen_14': 'Nitrate_reductionVIII(dissimilatory)',
						'nitrogen_15': 'Nitrate_reductionVIIIb(dissimilatory)',
						'nitrogen_16': 'Nitrate_reductionIX(dissimilatory)',
						'nitrogen_17': 'Nitrate_reductionX(dissimilatoryperiplasmic)',
						'nitrogen_18': 'Nitrogen_fixationI(ferredoxin)',
						'nitrogen_19': 'Nitrogen_fixation_II_(flavodoxin)',
						'nitrogen_20': 'Urea_degradationII',
						'nitrogen_21': 'Caffeine_degradationV(bacteria_via_trimethylurate)',
						'nitrogen_22': '4-aminobutanoate_degradationI',
						'nitrogen_23': '4-aminobutanoate_degradationII',
						'nitrogen_24': '4-aminobutanoate_degradationV',
						'nitrogen_25': 'allantoin_degradationIV(anaerobic)',
						'nitrogen_26': 'Ammonia_monoxygenase_AmoABC(Marker_gene)',
						'nitrogen_27': 'nirBD',
						'nitrogen_28': 'GABA_biosynthesis_prokaryotes_putrescine',
						'iron_1': 'Fe(II)oxidation',
						'iron_2': 'Fe_reduction_absorption',
                        'markers_1':'Single copy archaea',
                        'markers_2':'Single copy bacteria',
                        'markers_3':'Single copy both',
                        'markers_4':'miComplete bacteria',
                        'markers_5':'miComplete archaea',
                        }, inplace=True)

# Create a file with the genomic  completeness values
df_comp.to_csv(outpath/(args.filename + "_pfam_completenes.tab"), sep="\t")
df_gencomp=df_comp[['Single copy archaea','Single copy bacteria','Single copy both','miComplete bacteria','miComplete archaea']]
df_gencomp.to_csv(outpath/(args.filename + "_genomic_completenes.tab"), sep="\t")


# outfilename_comp = 'itol_mebs_comp.txt'
outfilename_comp = outpath/(args.filename + "_itol_mebs_comp.txt")
infile = "data2vis/dataset_heatmap_template.txt"
outfile2 = open(outfilename_comp, 'w')

# Modify FIELD_LABELS
with open(infile) as inf:
	for line in inf:
		if 'FIELD_LABELS ' in line:
			new_line = ['FIELD_LABELS'] + list(df_comp.columns)
			new_line = ' '.join(new_line) + '\n'   # spaces
			outfile2.write(new_line)
		else:
			outfile2.write(line)
# new line
outfile2.write('\n')

# new df

for ind_ in df_comp.index:
	l = [str(i) for i in df_comp.loc[ind_]]
	line = ind_ + ' ' + ' '.join(l) + '\n'
	outfile2.write(line)
outfile2.close()


# PLOTS
# Heatmap figure
sns.set(font_scale=0.7)
axs = sns.clustermap(df_comp.T, col_cluster=True, linewidths=0.1,
					 # cmap=sns.light_palette((210, 90, 60), input="husl"),
					 cmap=sns.diverging_palette(220, 20, n=10),
					 figsize=(15, 12))
# plt.tight_layout()
# plt.title("Metabolic completeness of S and C pathways", )
plt.savefig(outpath/(args.filename + "_comp_heatmap." + args.im_format),
			dpi=args.im_res, bbox_inches='tight')

# Barplot figure
if len(df.index) >= 50:
	print('''[WARNING] Skiping bar plot due to large data set (number of samples)''')
else:
	plt.figure(figsize=(12, 7))
	ax1 = df_new.plot(kind='bar')
	plt.ylabel("MEBS", weight='bold')
	plt.xlabel("Samples", weight="bold")
	plt.legend(bbox_to_anchor=(1, 1), loc=2,
			   borderaxespad=0.05, labelspacing=0.25)
	lines, labels = ax1.get_legend_handles_labels()
	ax1.legend(lines[:20], labels[:12], bbox_to_anchor=(
		1, 1), loc=2, borderaxespad=0.05, labelspacing=0.3)

	plt.savefig(outpath/(args.filename + "_barplot." + args.im_format),
				dpi=args.im_res, bbox_inches='tight')


# Heatmap
plt.figure(figsize=(7, 5))
ax1 = sns.heatmap(df_new.T,
				  #cmap=sns.color_palette("cubehelix", n_colors=10),
				  # cmap=sns.light_palette("dodgerblue", n_colors=10),
				  cmap=sns.diverging_palette(220, 20, n=10),
				  vmin=0, vmax=1, linewidths=0.1)
plt.ylabel("MEBS", weight='bold')
plt.xlabel("Samples", weight="bold")
plt.legend(bbox_to_anchor=(1, 1), loc=2, borderaxespad=0.05, labelspacing=0.25)
lines, labels = ax1.get_legend_handles_labels()
ax1.legend(lines[:20], labels[:12], bbox_to_anchor=(
	1, 1), loc=2, borderaxespad=0.05, labelspacing=0.3)
plt.savefig(outpath/(args.filename + "_mebs_heatmap." + args.im_format),
			dpi=args.im_res, bbox_inches='tight')

# Pairplot
df_new_ = df_new.copy()
df_new_ = df_new_.reset_index()
g = sns.PairGrid(df_new_.sort_values(["S", "C", "O", "Fe", "N"],
									 ascending=[False] * 5),
				 x_vars=["S", "C", "O", "Fe", "N"],
				 y_vars=["index"],
				 height=10, aspect=.25)
# Draw a dot plot using the stripplot function
g.map(sns.stripplot, size=10, orient="h",
	  palette="ch:s=1,r=-.1,h=1_r", linewidth=1, edgecolor="w")
# Use the same x axis limits on all columns and add better labels
g.set(xlim=(0, 1), xlabel="MEBS score", ylabel="")
# Use semantically meaningful titles for the columns
titles = ["S", "C", "O", "Fe", "N"]
for ax, title in zip(g.axes.flat, titles):

	# Set a different title for each axes
	ax.set(title=title)

	# Make the grid horizontal instead of vertical
	ax.xaxis.grid(False)
	ax.yaxis.grid(True)


sns.despine(left=True, bottom=True)
plt.tight_layout()
plt.savefig(outpath/(args.filename + "_mebs_dotplot." + args.im_format),
			dpi=args.im_res, bbox_inches='tight')

print("[END]   Done........................\n"
	  "Please check the following files in folder '{}'  :\n".format(str(outpath)),
	  "0. Original data without asteriks:", noast_fname, '\n',
	  "1. Heatmap displaying the metabolic completeness of N,Fe,S and CH4 pathways based on pfams:",
	  args.filename + "comp_heatmap.png\n",
	  "2. Barplot with normalized MEBS score values:", args.filename +
	  "_barplot.png\n",
	  "3. Heatmap with normalized MEBS score values:", args.filename +
	  "_mebs_heatmap.png\n",
	  "4. Dotplot with normalized MEBS score values:",  args.filename +
	  "_mebs_dotplot.png\n",
	  "5. Completeness file with description of the columns:", args.filename +
	  "_completenes.tab\n",
	  "6. Mapping file to itol with normalized MEBS scores:", args.filename + "_itol_mebs.txt\n",
	  "7. Mapping file to itol with pfam metabolic completeness:",  args.filename +
	  "_itol_mebs_comp.txt\n",
	  "8. File to be used as the input of F_MEBS_cluster.py -s none option", args.filename +
	  "_2_cluster_mebs.txt\n",
      "9. Genomic completeness based on single copy marker genes:", args.filename +
      " _genomic_completenes.tab\n",
      "10. Normalized mebs scores ", args.filename +
       "_norm_mebs.tab\n",
	  " If you have a tree file loaded in  itol, you can drag directly the _itol.txt files into your tree\n",
	  "and customize the colors of the pathways and the scores as in the following example\n",
	  "https://itol.embl.de/tree/97981518041461538630153\n",
      "Feel free to email me if you have any questions valdeanda[at]utexas[dot]edu"
	  )
