![MEBS](./images/MEBS.png) 

# Welcome to the new location of MEBS 

The development of  https://github.com/eead-csic-compbio/metagenome_Pfam_score   and any updates and new versions of MEBS will be posted on this new website.


# `Basic usage`

MEBS uses a few  command line options that can  be viewed by typing mebs.pl -h on the command line

```
perl mebs.pl -h 

  Program to compute MEBS for a set of genomic/metagenomic FASTA files in input folder.
  Version: v1.0

  usage: mebs.pl [options] 

   -help    Brief help message
   
   -input   Folder containing FASTA peptide files (.faa)                  (required)

   -type    Nature of input sequences, either 'genomic' or 'metagenomic'  (required)

   -fdr     Score cycles with False Discovery Rate 0.1 0.01 0.001 0.0001  (optional, default=0.01)

   -cycles  Show currently supported biogeochemical cycles
   
```

# About MEBS


 Multigenomic Entropy-Based Score or MEBS, allows the user to synthesizes genomic information into a single informative value. This entropy score can be used to infer the likelihood that microbial taxa perform specific metabolic-biogeochemical pathways

---

# `How does it work?`

<img src="https://valdeanda.github.io/mebs/images/mebs_overview.png"  align="right">

How can we organize and simplify complex genomic data to better understand the metabolisms of diverse microbial taxa? This was the major research question behind the development of MEBS. It works by synthesizing multiple sources of data on the metabolism of interest (e.g., genes, microbial taxa, and elemental reactions) under the mathematical framework of the Kullback-Leibler divergence, also known as relative entropy H’ (Fig. 1B). An H’ value is assigned to each protein-coding gene (pcg) based on how often the gene is encoded in the genomes of metabolically similar microbial taxa. H’ values close to 0 indicate that a given pcg is widely distributed among microorganisms (e.g., ABC transporters) and is therefore non-informative of the target metabolism. H’ values near or greater than 1 indicate that a given pcg is unique to a group of metabolically similar microbial taxa, whereas negative H’ values indicate a given pcg is not expected to be involved in the target metabolism. MEBS stores the unique pcgs and their H’ values in an internal database, which can then be cross-referenced with genomic and metagenomic input data to obtain a single H’ score for a given metabolic pathway. High H’ scores suggest that microbial taxa can perform the pathways involved in the metabolism of interest.

Currently, the MEBS software has a built-in function that allows users to accurately and quickly evaluate the likelihood that up to thousands of microbial taxa (genomes) or communities (metagenomes) can perform the metabolic reactions involved in C, N, O, S, and Fe cycling. However, a user can also take a more advanced, step-wise, approach to examine additional metabolic reactions (e.g., As cycling).

In addition, MEBS includes several tools to help users interpret their results: H’ values of specific metabolic pathways can be mapped on visual reconstructions of microbial phylogenies (Fig. 1C); microbial taxa or entire communities can be sorted by their metabolic function (e.g., nutrient assimilation) (Fig. 1D); or communities filling certain metabolic niches (e.g., heavy metal degradation) can be mapped at global scales (Fig. 1E).

---
 
# `Aplications` 

MEBS provides an open-access, reproducible, entropy-based platform to efficiently analyze microbial genomic datasets and decode a plethora of metabolic-chemical pathways. The broad application of MEBS stretches across scientific disciplines. MEBS can be used to study microbial responses to on-going threats of anthropogenic and climate-linked change associated with sea level rise, ocean acidification, and land-use by providing new information on shifts in microbial-mediated biogeochemical cycling and bottom-up effects on entire ecosystems. MEBS can identify global environments with microbial taxa that carry the necessary pcgs for the degradation of hydrocarbons, heavy metals, and other pollutants, which lends itself to new bioremediation practices. MEBS can also be used in translational medicine to develop new tests that differentiate groups of infectious, antibiotic-resistant, bacteria in the human lung. In the future, I see MEBS software along with other metabolic decoder approaches integrated into sequencer instruments (e.g., MinION), allowing real-time analysis of complex metabolic pathways in nature


# `Real life applications of MEBS`

## Protein clustering using mebs  
<img src="https://valdeanda.github.io/mebs/images/deltas.png" width="200" height="150" align="right">

 Group 2K metabolically-related Deltaproteobacteria genomes using mebs_clust.py script 
 See [Langwig-De Anda et al. 2021](https://www.nature.com/articles/s41396-021-01057-y)










# Dependencies

- python >= 3.6.2
- python modules:
  - numpy
  - matplotlib
  - pandas
  - seaborn

















































