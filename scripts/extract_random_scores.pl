#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin '$Bin';

# This script takes a folder and reads all score files therein with the following pattern:
#4489656.3.27122016.hmmsearch.tab.17.score
#       <- hmmsearch id 
#               <- output date 
#                                  <-replicate
#or
#4517710.3.hmmsearch.tab.99.score
#       <- hmmsearch id 
#                        <-replicate
#or 
#GCF_000211375.1_ASM21137v1_protein.faa.named.faa.out.hmmsearch.tab.163.score
#
#   <-genome    

# Output: 
# 1) TAB-separated output with parsed scores 
# #columns replicates 
# #rows genomes or metagenomes 
# B Contreras-Moreira, V de Anda 2016

my $DEFAULTDIR = $Bin.'/scores';
my $FILENAMEPATTERN = '(\S*?\d+\.\d+).*?\.hmmsearch\.tab\.(\d+).score';
# aqui hay varios tipos de archivos todos
# lo mejor es que el nombre para los genomas sean los primeros 14 digitos y para metagenomas los primeros 8 
# la replica siempre será la misma *(\d+).score  
#
#
my ($INP_help,$INP_dir) = ('',$DEFAULTDIR);

GetOptions
(
    'help|h|?'        => \$INP_help,
	 'dire|dir=s' => \$INP_dir,
);

if (-t STDIN && $INP_help )
{
die<<EODOC;

Program to merge Pfam entropy scores from multiple random experiments, which can then be plotted.

usage: $0 [options] 

 -help             brief help message
 
 -dire        directory containing Pfam entropy scores from random experiments  (string, \n   default $DEFAULTDIR)
 
EODOC
}

$INP_dir =~ s/\/$//;
   
print "# $0 call:\n# -matrixdir $INP_dir\n\n";

################################################

my ($id,$replicate,$ss,%scores);
my $max_replicate = 0;

## 1) find matrix files

opendir(DIR,$INP_dir) || die "# $0 : ERROR : cannot list $INP_dir\n";
my @files = grep{/$FILENAMEPATTERN/} readdir(DIR);
closedir(DIR);

## 2) parse all files
foreach my $file (@files)
{
  open(INFILE,"$INP_dir/$file") || die "# $0 : cannot find $INP_dir/$file\n";

  # extract name and replicate
  if($file =~ m/$FILENAMEPATTERN/)
  {
    ($id,$replicate) = ($1,$2);
  }

  while(<INFILE>)
  { 	
    #Pfam entropy score: 3.645
    #Pfam entropy score: 0
    if(/Pfam entropy score:\s+(\S+)/)
    {
      $scores{$id}[$replicate] = $1;      
      if($replicate > $max_replicate){ $max_replicate = $replicate }
    }
  }
  close(INFILE);
} 
	

## produce merged score files 

open(OUTFILE,'>',$INP_dir.".random_scores.tab");

# header
foreach $replicate (1 .. $max_replicate)
{ 
  print OUTFILE "\t$replicate";
}  print OUTFILE "\n";

foreach $id (sort keys(%scores))
{
  print OUTFILE "$id";

  foreach $replicate (1 .. $max_replicate) 
  {
    if(defined($scores{$id}[$replicate])){ $ss = $scores{$id}[$replicate] }
    else { $ss = 'NA' }
    print OUTFILE "\t$ss";
  } 
  print OUTFILE "\n";
  
}

close(OUTFILE);





