![MEBS](./images/MEBS.png) 

# Welcome to the new location of MEBS 

The development of  https://github.com/eead-csic-compbio/metagenome_Pfam_score   and any updates and new versions of MEBS will be posted on this new website.


# Dependencies

- python >= 3.6.2
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
## Protein clustering using mebs  

 Group 2K metabolically-related Deltaproteobacteria genomes using   [mebs_clust.py](https://github.com/valdeanda/mebs/blob/master/mebs_clust.py) [Langwig-De Anda et al. 2021](https://www.nature.com/articles/s41396-021-01057-y)

```
 python3 mebs_clust.py  -h
usage: mebs_clust.py [-h] [--outdir OUTDIR]
                     [-im_format {png,pdf,ps,eps,svg,tif,jpg}] [--im_res dpi]
                     [--cutoff CUTOFF]
                     [--method {single,complete,average,weighted,centroid,median,ward}]
                     [--distance {braycurtis,canberra,chebyshev,cityblock,correlation,euclidean,jaccard,mahalanobis}]
                     [--nolegend]
                     filename

positional arguments:
  filename              Input file pfam profile.

optional arguments:
  -h, --help            show this help message and exit
  --outdir OUTDIR, -o OUTDIR
                        Output folder [<filename>_mebs_clust]
  -im_format {png,pdf,ps,eps,svg,tif,jpg}, -f {png,pdf,ps,eps,svg,tif,jpg}
                        Output format for images [pdf].
  --im_res dpi, -r dpi  Output resolution for images in dot per inch (dpi)
                        [dpi].
  --cutoff CUTOFF, -co CUTOFF
                        Cutoff threashold default 0.5
  --method {single,complete,average,weighted,centroid,median,ward}, -m {single,complete,average,weighted,centroid,median,ward}
                        methods for calculating the distance between the newly
                        formed clusters [ward].
  --distance {braycurtis,canberra,chebyshev,cityblock,correlation,euclidean,jaccard,mahalanobis}, -d {braycurtis,canberra,chebyshev,cityblock,correlation,euclidean,jaccard,mahalanobis}
                        The distance metric to use default [jaccard].
  --nolegend, -nl       If specified, the legend is not shown

Example:
$  python3 mebs_clust.py  pfamprofile.tsv
``` 


<img src="https://valdeanda.github.io/mebs/images/deltas.png" width="200" height="150" align="right">





























































