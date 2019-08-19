#!/usr/bin/perl -w
use strict;
use Getopt::Long; 
#updated 2018 
 ## This script takes 1-2 inputs:
 ## 1) a directory containing the  hmmsearch TSV outfiles with the results of scanning a collection metagenomes against
 ## Pfam-A database
 ## Output:
 ## 1) a matrix of occurrence of Pfam domains across the metagenomic sample/ raw/ normalized / binary ? 
 #
 ## B Contreras-Moreira, V de Anda 2017
 
 
my $COUNTPFAMONCEPERSEQ = 1; # if set to zero counts all, repeated instances of Pfam domains per sequences

my $FILENAMEPATTERN = 'S(\d+)\.pfam.tab';

my ($INP_help,$INP_matrixdir) = ('','');

GetOptions(
  'help|h|?'        => \$INP_help,
  'matrixdir|dir=s' => \$INP_matrixdir,
);

if (-t STDIN && $INP_help )
{
die<<EODOC;
Program to merge hmmsearch tab-separed outfiles into a single matrix of Pfam ocurrence across multiple samples

usage: $0 [options] 

 -help         brief help message
 
 -matrixdir    directory containing tsv output files from hmmsearch ze (string) (required)

 -pfam_list    tabular format file containing the optional Pfams to compute the ocurrence matrix (optional) 

EODOC
}

if (!$INP_matrixdir)
{
die "#Error: require valid -matrixdir directory. Type -h to see the options \n";
}



$INP_matrixdir =~ s/\/$//;
   
print "# $0 call:\n# -matrixdir $INP_matrixdir\n";
print "# PARAMS: COUNTPFAMONCEPERSEQ=$COUNTPFAMONCEPERSEQ\n\n";

################################################
#


# 1) find matrix files
opendir(MATDIR,$INP_matrixdir) || die "# $0 : ERROR : cannot list $INP_matrixdir\n";
my @matfiles = grep{/$FILENAMEPATTERN/} readdir(MATDIR);
closedir(MATDIR);

# 2) parse all files
my ($line,$sample_id,$seqid,$domid,%hmm,@samples);
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

  #3)  extract sample id   from filename
  if($matrixfile =~ m/S(\d+).pfam.tab/)
  {
    $sample_id = $1;
    push(@samples,$sample_id);
  }
  else{ die "# cannot retrieve sample_id from $matrixfile\n"; }


  #5) read hmmsearch output file (tblout), last column in [] should be the the input name
  my %countsperseq;
  while($line = <INFILE>)
  {
    # sample contents:
    # target name accession query name accession ...  description of target
    #WP_048201923.1       -          2-Hacid_dh           PF00389.28  ... [GCF_000739065.1_ASM73906v1_protein.faa]
    next if($line =~ /^#/); #print $line;

    my @data = split(/\s+/,$line);

    ($seqid,$domid) = ($data[0],$data[3]);

    $domid = (split(/\./,$domid))[0];

    $countsperseq{$seqid}{$domid}++;

    next if($COUNTPFAMONCEPERSEQ && 
	$countsperseq{$seqid}{$domid} && $countsperseq{$seqid}{$domid} > 1);

    $hmm{$domid}{$sample_id}++; #print "$sample_id $domid $hmm{$domid}{$sample_id}\n";
  }
  close(INFILE);
}

#6) Print the count matrix
my @sorted_domains = sort(keys(%hmm));

foreach $sample_id (@samples)
{
	print "\t$sample_id";
} print "\n";
	
foreach $domid (@sorted_domains)
{
	print "$domid";
	foreach $sample_id (@samples)
	{
		print "\t$hmm{$domid}{$sample_id}";
	} print "\n";
}




