#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin '$Bin';


# This script takes 2-3 inputs:
# 1) a hmmsearch TSV outfile with the results of scanning a collection of Pfam domais against
# a large set of (non-redundant) genomes
# 2) a list of selected accessions of genomes interest to compute entropy
# 3) an optional list of RefSeq assembly annotations to print scientific names instead of accession codes
#
# Output:
# 1) a matrix of occurrence of Pfam domains across genomes
# 2) entropy estimates of each scanned Pfam domain with respect to the selected accessions
# B Contreras-Moreira, V de Anda 2016
# Last updated August 2018

my $INP_help = 0;
my $INP_absolute_freqs = 0;
my $PSEUDOCOUNT = 0.8; #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2647310/
my ($dom_assign_file,$list_file,$refseq_file) = ('','','');


GetOptions
(
  'help|h|?'        => \$INP_help,
  'input_dom|in=s'  => \$dom_assign_file,
  'input_list|l=s'  => \$list_file,
  'names|n=s'       => \$refseq_file,
  'absolute|a'      => \$INP_absolute_freqs
);

if (-t STDIN && ($INP_help || $dom_assign_file  eq '' || $list_file eq ''))
{
die<<EODOC;

Program to compute Pfam entropies from a list of accession genomes of interest.

usage: $0 [options] 

 -help          Brief help message
 
 -input_dom     tbl-format file with HMM matches produced by hmmsearch (required)

 -input_list    list of selected genomes of interest                   (required) 

 -names         optional list of Refseq assembly annotations to print
                scientific names instead of accesion codes             (optional) 

 -absolute      print absolute domain frequencies                      (optional, default binary)

EODOC
}




#Required input files 
if(!$dom_assign_file || !-s $dom_assign_file)
{ 
 die "#ERROR: cannot locate input file -input_dom $dom_assign_file\n";

}
elsif (!$list_file || !-s $list_file)
{
 die "ERROR: cannot locate input file -input_list $list_file\n";
}

#Optional arguments 

if ($refseq_file && !-s $refseq_file)
{
 die "ERROR: cannot locate -names file $refseq_file\n";
}


print "# $0 call:\n# -input_dom $dom_assign_file -input_list $list_file  -names $refseq_file\n\n";

#################################

#if(!$ARGV[1]){ die "# $0 : usage: $0 <pfam_hmmsearch.tab> <accession list (ie Suli)> <RefSeq list, optional>\n"; }
#else
#{ 
#  ($dom_assign_file,$list_file,$refseq_file) = @ARGV;
#  if(!$refseq_file){ $refseq_file = 'null' } 
#}

#print "# dom_assign_file=$dom_assign_file\n# list_file=$list_file\n# RefSeq_file=$refseq_file\n";
print "# PSEUDOCOUNT=$PSEUDOCOUNT\n\n";




my ($seqid,$domid,$taxonid,$spid,$line,$total_obs,$total_exp,$freq_obs,$freq_exp,$entropy);
my (%taxa,%hmm,%domain_seqid,%matrix,%list_matrix,%taxon_list,%matched_taxa,%scnames);
my $listOK = 0;

## 0) read optional RefSeq file
if($refseq_file && -s $refseq_file)
{
  open(REFSEQ,$refseq_file) || die "# $0 : cannot read $refseq_file\n";
  while($line = <REFSEQ>)
  {
    #GCF_000262325.2  PRJNA224116 SAMN02604280    na  1037911 294 Pseudomonas fluorescens A506  strain=A506   latest...
    next if($line =~ /^#/); 
    my @data = split(/\t/,$line);
    $scnames{$data[0]} = $data[7]; #print "$data[0] $data[7]\n";
  }
  close(REFSEQ);
}

## 1) read accession list file
if($list_file && -s $list_file)
{
    open(LIST,$list_file) || die "# $0 : cannot read $list_file\n";
    while(<LIST>)
    {
        #List of accessions & genomes of interest
        #GCF_000213215.1 Acidianus hospitalis W1
        #GCF_000020825.1 Acidithiobacillus ferrooxidans ATCC 53993
        #GCF_000175575.2 Acidithiobacillus caldus ATCC 51756

        chomp;
        $taxonid = (split(/\s+/,$_))[0];
        #next if ($taxonid =~ /^#/);
        $taxonid =~ s/\s+$//g;
        $taxonid =~ s/\.$//g; #next if($taxonid =~ /sp\./ || /^[a-z]/);

        $taxon_list{$taxonid}++;
    }
    close(LIST);

    printf("# number of taxa in %s = %d\n\n",$list_file,scalar(keys(%taxon_list)));
    $listOK = 1;
}

## 2) read hmmsearch output file (tblout), last column in [] should be the the input name
open(DOMFILE,$dom_assign_file) || die "# $0 : cannot read $dom_assign_file\n";
while($line = <DOMFILE>)
{
    # sample contents:
    # target name accession query name accession ...  description of target
    #WP_048201923.1       -          2-Hacid_dh           PF00389.28  ... [GCF_000739065.1_ASM73906v1_protein.faa]
    next if($line =~ /^#/); #print $line;

    my @data = split(/\s+/,$line);

    ($seqid,$domid) = ($data[0],$data[3]);

    $domid = (split(/\./,$domid))[0];

    $domain_seqid{$seqid}{$domid}++;

    # count domains once per sequence
    next if($domain_seqid{$seqid}{$domid} > 1);

    if(!$hmm{$domid}){ $hmm{$domid}++ }

    # find a bracketed [accession name/number] 
    $line = reverse($line);
    if($line =~ /\](.+?)\[/)
    {
        $taxonid = reverse($1);
        if($taxonid =~ m/(\w+_\d+\.\d+)/){ $taxonid = $1 } # such as GCF_000520015.2

        # adhoc rules to filter out common annotations such as:
        # [5'-phosphate] , [3-hydroxymyristoyl] ,
        # [acyl-carrier-protein] , [NAD(P)H] , [Cu-Zn]
        if(!$taxa{$taxonid}){ $taxa{$taxonid}++ }
    }
    else{ next; }

    # check wether this taxonid is in $list_file
    if($listOK && $taxon_list{$taxonid})
    {
      $list_matrix{$taxonid}{$domid}++;    
    }

    # save the occurrence of this domain/HMM in this taxon
    $matrix{$taxonid}{$domid}++; 
}
close(DOMFILE);

## 3) Print MATRIX with presence/absence of PFAM domains

# Print Header
foreach $domid (sort keys(%hmm))
{
    print "\t$domid";
} print "\n";

#Print Data 
foreach $taxonid (sort keys(%taxa))
{
    printf("%s",$scnames{$taxonid} || $taxonid);

    foreach $domid (sort keys(%hmm))
    {
        if($matrix{$taxonid}{$domid})
        {
          if($INP_absolute_freqs){ print "\t$matrix{$taxonid}{$domid}" }
          else{ print "\t1" }

          $matrix{'total'}{$domid}++;
          if($list_matrix{$taxonid})
          {
            $list_matrix{'total'}{$domid}++;
            $matched_taxa{$taxonid}++;
          }
        }
        else{ print "\t0" }
    } print "\n";
}

#Print expected frequency for each domain  (all dataset ) )
printf("\n\nexpfreq(%d)",scalar(keys(%taxa))-1);
foreach $domid (sort keys(%hmm))
{
    $total_exp = $matrix{'total'}{$domid} || 0; 
    printf("\t%1.3f",$total_exp/(scalar(keys(%taxa))-1));
} print "\n";

## 4) Compute the relative entropy 

#Observed frequency P(i): occurrences of protein family i in sulfur-related genomes
if($listOK)
{
    # Print the organisms matched in the list
    printf("\n\nmatched list taxa in %s (%d) : ",
      $dom_assign_file,
      scalar(keys(%matched_taxa)));

    foreach $taxonid (sort keys(%matched_taxa))
    {
        printf("%s,",$scnames{$taxonid} || $taxonid);
    }
    print "\n";

    # Print P(i) for each Pfam in the input list 
    printf("\n\nlistfreq(%d)",scalar(keys(%matched_taxa)));
    foreach $domid (sort keys(%hmm))
    {
        $total_obs = $list_matrix{'total'}{$domid} || 0;
        printf("\t%1.3f",$total_obs/scalar(keys(%matched_taxa)));
    }
    print "\n";

    #Print relative entropies, in bits
    #Entropy captures the extent to which a family informs specifically about specific metabolisms (such as S metabolism)
    #Largest values correspond to the most informative families (sulfur-energy based genomes, whereas values close to zero
    #describe non-informative families
    #Negative values correspond to protein families observed less than expected, and might be also markers
    printf("\n\nrel_entropy(%d)",scalar(keys(%matched_taxa)));
    foreach $domid (sort keys(%hmm))
    {   
        $total_exp = $matrix{'total'}{$domid} || $PSEUDOCOUNT;
        $freq_exp = $total_exp/(scalar(keys(%taxa))-1);

        $total_obs = $list_matrix{'total'}{$domid} || $PSEUDOCOUNT; 
        $freq_obs = $total_obs/scalar(keys(%matched_taxa));

        if($freq_obs && $freq_exp)
        {
            $entropy = sprintf("%1.3f",$freq_obs * (log($freq_obs/$freq_exp)/log(2)));
        }
        else{  $entropy = 'NA' }

        print "\t$entropy";
    }   
    print "\n";
}
