![MEBS](./images/MEBS.png) 

# Welcome to the new location of MEBS 

The development of  https://github.com/eead-csic-compbio/metagenome_Pfam_score   and any updates and new versions of MEBS will be posted on this new website.

# Dependencies


- hmmsearch 
- python >= 3.6.2
- python modules:
  - numpy
  - matplotlib
  - pandas
  - seaborn


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
The code is regularly patched (see [CHANGES.txt](./CHANGES.txt)regurarly. 

MEBS has been used in a variety of studies (see citing papers [here](https://scholar.google.com/scholar?oi=bibs&hl=en&cites=4642441397530015315)

We kindly ask you to report errors or bugs in the program to the authors and to acknowledge the use of the program in scientific publications.


 Installation instructions are summarized on [README.txt](./README.txt) and full documentation is available in two        flavours:

 |version|Link|
 |-------|----|
 |Updated MEBS documentation|[manual](https://valdeanda.github.io/mebs/README.html)|
 |First MEBS documentation in pdf|[manual-pdf](https://github.com/eead-csic-compbio/metagenome_Pfam_score/blob/master/manual.v1.pdf/)|

