#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin '$Bin';

# This script takes a folder and reads all entropy files therein with the following pattern:
# random1.txt.100.csv 
#       <- number of random iteration
#               <- peptide size
#
# Infiles might be compressed with GZIP or BZIP2
#
# Output: 
# 1) TAB-separated output with Pfam entropies of parsed files

# B Contreras-Moreira, V de Anda 2016

my $DEFAULTMATRIXDIR = $Bin.'/../data/entropies_matrix/';
my $FILENAMEPATTERN = 'random(\d+)\.txt.*?csv';

my ($INP_help,$INP_matrixdir) = ('',$DEFAULTMATRIXDIR);

GetOptions
(
    'help|h|?'        => \$INP_help,
	 'matrixdir|dir=s' => \$INP_matrixdir,
);

if (-t STDIN && $INP_help )
{
die<<EODOC;

Program to merge Pfam entropy scores from multiple experiments, which can then be plotted.

usage: $0 [options] 

 -help             brief help message
 
 -matrixdir        directory containing pre-computed random entropies from peptides of variable size (string, \n                    default $DEFAULTMATRIXDIR)
 
EODOC
}

$INP_matrixdir =~ s/\/$//;
   
print "# $0 call:\n# -matrixdir $INP_matrixdir\n\n";

################################################

my (@HMMentropy,@HMMs,%sizes);
my ($hmm,$entropy,$replicate,$size);

## 1) find matrix files
opendir(MATDIR,$INP_matrixdir) || die "# $0 : ERROR : cannot list $INP_matrixdir\n";
my @matfiles = grep{/$FILENAMEPATTERN/} readdir(MATDIR);
closedir(MATDIR);

## 2) parse all files
foreach my $matrixfile (@matfiles)
{
  if($matrixfile =~ m/\.bz2/)
  {
    open(INFILE,"bzcat $INP_matrixdir/$matrixfile|") || die "# $0 : cannot find $INP_matrixdir/$matrixfile\n";
  }
  elsif($matrixfile =~ m/\.gz/)
  {
    open(INFILE,"zcat $INP_matrixdir/$matrixfile|") || die "# $0 : cannot find $INP_matrixdir/$matrixfile\n";
  }
  else
  {
    open(INFILE,"$INP_matrixdir/$matrixfile") || die "# $0 : cannot find $INP_matrixdir/$matrixfile\n";
  }

  # extract peptide size from filename
  if($matrixfile =~ m/random(\d+)\.txt\.(\d+)\.csv/)
  {
    ($replicate,$size) = ($1,$2);
  }
  elsif($matrixfile =~ m/random(\d+)\.txt/)
  {
     $replicate = $1;
     $size = 'real';
  }
  
  # store sizes observed
  $sizes{$size}++;
  
      
  while(<INFILE>)
  { 	
  	#	PF00005	PF00009	...
	  if(/^\tPF\d+\t/)
	  {
		  chomp;
		  @HMMs = split(/\t/,$_); 
      shift(@HMMs); # delete first cell, empty
	  }
	  elsif(/^rel_entropy/)
	  {
		  chomp;
		  my @entropies = split(/\t/,$_);
      shift(@entropies);
		  foreach $hmm (0 .. $#HMMs)
		  {
			  $HMMentropy[$replicate]{$HMMs[$hmm]}{$size} = $entropies[$hmm];
		  }
		
		  last;
	  }
  }
  close(INFILE);
}
	

## produce merged score files 
foreach $size (sort keys(%sizes))
{
  open(OUTFILE,'>',$INP_matrixdir.".$size.tab");

  print "# merged file of MSL=$size (replicates=$sizes{$size}, $INP_matrixdir.$size.tab)\n";

  # header
  foreach $replicate (1 .. $#HMMentropy)
  {
    print OUTFILE "\t$replicate";
  }  print OUTFILE "\n";

  foreach $hmm (0 .. $#HMMs)
  {
    print OUTFILE "$HMMs[$hmm]";

    foreach $replicate (0 .. $#HMMentropy)
    {
      $entropy = $HMMentropy[$replicate]{$HMMs[$hmm]}{$size} || 'NA';
      print OUTFILE "\t$entropy";
    } print OUTFILE "\n";
  }

  close(OUTFILE);
}  




