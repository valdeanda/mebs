![MEBS](./images/MEBS.png) 


# Dependencies

- python >= 3.4
- python modules:
  - numpy
  - matplotlib
  - pandas
  - seaborn


# Basic usage

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































































