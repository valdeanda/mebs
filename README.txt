Authors
===============================================================================
MEBS  is designed, created and maintained by Valerie De Anda at the University of Austin Texas, USA (https://cns.utexas.edu/component/cobalt/item/9-marine-science/3799-de-anda-valerie?Itemid=349) Bruno Contrearas Moreira at the Laboratory of Computational
Biology at Estacion Experimental de Aula Dei (EEAD/CSIC http://www.eead.csic.es/compbio)
in Zaragoza, Spain, and Cesar Augusto Hernandez at the Cellular Physiology Institute of Universidad Nacional
Autonoma de Mexico (http://www.ifc.unam.mx/busqueda.php).

The program was written mostly by Valerie De Anda, Bruno Contreras-Moreira and Cesar Augusto Hernandez,
but includes code and binaries from other authors:

 Perl5  (https://www.perl.org/get.html), preinstalled in most Linux systems
 HMMER 3.1b2 (http://hmmer.org)
 Pfam (http://pfam.xfam.org , PubMed=19920124)
Key literature references:
================================================================================

1. De Anda V, Zapata-Peñasco I, Poot-Hernandez AC, Eguiarte LE, Contreras-Moreira B, Souza V. MEBS, a software platform to evaluate large (meta)genomic collections according to their metabolic machinery: unraveling the sulfur cycle. Gigascience. 2017 Nov 1;6(11):1-17. doi: 10.1093/gigascience/gix096. PMID: 29069412; PMCID: PMC5737871.

2.De Anda V, Zapata-Peñasco I, Blaz J, Poot-Hernández AC, Contreras-Moreira B, González-Laffitte M, Gámez-Tamariz N, Hernández-Rosales M, Eguiarte LE, Souza V. Understanding the Mechanisms Behind the Response to Environmental Perturbation in Microbial Mats: A Metagenomic-Network Based Approach. Front Microbiol. 2018 Nov 28;9:2606. doi: 10.3389/fmicb.2018.02606. PMID: 30555424; PMCID: PMC6280815.

3.Langwig MV, De Anda V, Dombrowski N, Seitz KW, Rambo IM, Greening C, Teske AP, Baker BJ. Large-scale protein level comparison of Deltaproteobacteria reveals cohesive metabolic groups. ISME J. 2021 Jul 30. doi: 10.1038/s41396-021-01057-y. Epub ahead of print. PMID: 34331018.




System requirements and installation:
================================================================================
This software has been tested in Linux and MacOSX systems. The software is composed of 3 main scripts (mebs.pl, mebs_vis.py and mebs_clust.py) written in perl, python3 and python3 respectively.

1. mebspl : In order to use mebs.pl you only need to have perl and hmmsearch installed

To install hmmsearch see documentation webpage to obtain and compile HMMER from source:
http://hmmer.org/documentation.html

  $  brew install hmmer               # OS/X, HomeBrew
  $  port install hmmer               # OS/X, MacPorts
  $  apt install hmmer                # Linux (Ubuntu, Debian...)
  $  dnf install hmmer                # Linux (Fedora)
  $  yum install hmmer                # Linux (older Fedora)
  $  conda install -c bioconda hmmer  # Anaconda

2. mebs_vis.py: For the visualization script you will need  python >= 3.6.2 and several python modules
        numpy
        matplotlib
        pandas
        seaborn

To install the above modules including python3

   $ sudo apt-get install python3             # Linux (Ubuntu, Debian..)
   $ sudo pip3 install -U pip                 # Linux (Ubuntu, Debian..)
   $ sudo -H pip3 install --upgrade pandas    # Linux (Ubuntu, Debian..)
   $ sudo -H pip3 install --upgrade numpy     # Linux (Ubuntu, Debian..)
   $ sudo -H pip3 install --upgrade scipy     # Linux (Ubuntu, Debian..)
   $ sudo -H pip3 install --upgrade seaborn   # Linux (Ubuntu, Debian..)


Conda option

Istall conda or miniconda depending of your needs. See webpage for more info
https://conda.io/projects/conda/en/latest/user-guide/install/index.html

Once conda is installed activate conda and use the yml file provided with mebs.

conda activate
conda env create -f mebs.yml
conda activate mebs_env


