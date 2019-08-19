#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin '$Bin';

# This script takes two input files and several optional choices: 
#
# 1) a hmmsearch tbl-format outfile with the results of scanning a collection of Pfam domais against
# a large set of genomes or metagenomes
#
# 2) a TSV file with pre-computed entropies of Pfam domains produced by scripts/extract_entropies.py

# Output:
# 1) TAB-separated output with Pfam entropy scores

# B Contreras-Moreira, V de Anda 2017

my $DEFAULTFRAGSIZE  = 100;
my $DEFAULTMINRELENTROPY  = -9;

# RGB scale from pink to red (unmatched KEGG nodes should be white #FFFFFF)
#HEX scale from blue, yellow to red 
my @COLORS = ( '#2c7bb6', '#abd9e9',  '#ffffbf', '#fdae61', '#d7191c' );
#my @COLORS = ( '#FF9999', '#FF6666',  '#FF3333', '#FF0000', '#CC0000' );

my ($INP_pfamsearchfile,$INP_infile_bzipped,$INP_entropyfile) = ('',0,'');
my ($INP_fragment_size,$INP_help,$INP_keggmapfile,$INP_minentropy) = ($DEFAULTFRAGSIZE,0,'',$DEFAULTMINRELENTROPY);
my ($INP_pathways,$RAND_percent,@user_pathways,$pw) = ('',0);

GetOptions
(
  'help|h|?'        => \$INP_help,
  'input|in=s'      => \$INP_pfamsearchfile,
  'size|s:s'        => \$INP_fragment_size,
	'bzip|b'          => \$INP_infile_bzipped,
	'entropyfile|enf=s' => \$INP_entropyfile,
  'keggmap|km=s'    => \$INP_keggmapfile,
  'minentropy|min=f'=> \$INP_minentropy,
	'pathway|pw=s'    => \$INP_pathways,
	'random|rd:i'     => \$RAND_percent
);

if (-t STDIN && ($INP_help || $INP_pfamsearchfile eq '' || $INP_entropyfile eq ''))
{
die<<EODOC;

Program to compute Pfam entropy score (i.e Sulfur score).

usage: $0 [options] 

 -help          Brief help message
 
 -input         tbl-format file with HMM matches produced by hmmsearch (required)

 -size          Mean peptide size length (MSL)                         (optional, default=$INP_fragment_size)
                MSL takes integers 30,60,100,150,200,250,300 for the 
                analysis of metagenomes, where peptides are usually 
                fragmented, and also the value 'real' in the case of
                completely sequenced genomes.

 -bzip          Input file is BZIP2-compressed                         (optional)
 
 -entropyfile   TSV file with pre-computed entropies from peptides     (required)
                of variable size, produced by extract_entropies.py
                The file must be formatted as follows:

                  real  30  60  100 150 200 250 300
                  PF00005 -0.001  0.001 -0.001  -0.001  -0.001  -0.001  -0.001  -0.001
                  ...
 
 -minentropy    Min relative entropy of HMMs to be considered          (optional, unsigned float, example -min 0.3)
                to compute the Pfam score    
						 
 -keggmap       TSV file with HMM to KEGG mappings and pathway names.  (optional)
                This produces a report on pathway completness and  
                also a script to display the pathways in KEGG.
                The file must be formatted as follows: 

                  PFAM  KO  PATHWAY   PATHWAY NAME 
                  PF00890 K00394  1 Sulfite oxidation 
                  PF01087   1 Sulfite oxidation 
                  PF00581 K01011  2 Thiosulfate oxidation
                  ...
 
 -pathway       Comma-separated pathway numbers from -keggmap file to  (optional, by default all pathways are used)
                produce pathway reports and to compute the Pfam score

 -random        Percent of random-sampled Pfams to compute the score   (integer, default=100)

EODOC
}
   
# required input files
if(!$INP_pfamsearchfile || !-s $INP_pfamsearchfile)
{
    die "# ERROR : cannot locate input file -input $INP_pfamsearchfile\n";
}  
elsif(!$INP_entropyfile || !-s $INP_entropyfile)
{
  die "# ERROR : cannot locate input file -entropyfile $INP_entropyfile\n";
}

# optional arguments
if($INP_fragment_size ne 'real' && ($INP_fragment_size !~ /^\d+$/ || $INP_fragment_size < 1))
{
    die "# ERROR : invalid value for fragment size ($INP_fragment_size)\n";
}

if($INP_keggmapfile && !-s $INP_keggmapfile)
{
    die "# ERROR : cannot locate input file -keggmap $INP_keggmapfile\n";
}

if($RAND_percent && ($RAND_percent < 1 || $RAND_percent > 100))
{
    die "# ERROR : invalid percent to compute the scores ($RAND_percent)\n";
}

elsif($INP_keggmapfile)
{
	if($INP_pathways)
	{
		foreach $pw (split(/,/,$INP_pathways))
		{
			push(@user_pathways,$pw);
		}
	}
}

print "# $0 call:\n# -input $INP_pfamsearchfile -size $INP_fragment_size -bzip $INP_infile_bzipped ".
	"-entropyfile $INP_entropyfile -minentropy $INP_minentropy -random $RAND_percent -keggmap $INP_keggmapfile ".
  "-pathway $INP_pathways\n\n";

################################################

my ($random,%skip_random,$r) = ( $RAND_percent ); 

my (%HMMentropy,@HMMs,%matchedHMMs,%KEGGmap,$hmm,$KEGGid,$entropy);
my $size_column = -1;

## open entropy file, validate the selected size, and store entropies
open(ENTROPYFILE,$INP_entropyfile) || die "# $0 : cannot find $INP_entropyfile\n";
while(<ENTROPYFILE>)
{
  chomp;

  # header file with accepted fragment sizes 
  if(/^\treal/) # real  30  60  100 150 200 250 300
  {
    my @valid_sizes = split(/\t/,$_);
    shift(@valid_sizes); # first column is empty

    foreach my $vsize (0 .. $#valid_sizes)
    {
      if($valid_sizes[$vsize] eq $INP_fragment_size)
      {
        $size_column = $vsize;
        last;
      }
    }

    if($size_column == -1)
    {
      die "# ERROR : cannot find -size $INP_fragment_size in -entropyfile $INP_entropyfile\n";
    }
  }
  else # PF00005 -0.001  0.001 -0.001  -0.001  -0.001  -0.001  -0.001  -0.001
  {
    my @entropies = split(/\t/,$_);

    # save ordered list of parsed HMMs 
    $hmm = shift(@entropies);
    push(@HMMs,$hmm);  

    # store selected entropy according to size
    $HMMentropy{$hmm} = $entropies[$size_column];
  }
}
close(ENTROPYFILE);

printf("# total HMMs with assigned entropy in %s : %d\n\n",
  "$INP_entropyfile",scalar(keys(%HMMentropy)));

## if requested, random-sample a fraction of parsed HMMs
if($RAND_percent) 
{ 
  my $hmm_number = scalar(@HMMs);
  my $total_hmms2skip = int( ((100-$random)/100) * $hmm_number );
        
  while(scalar(keys(%skip_random)) < $total_hmms2skip)
  {
    $r = int(rand($#HMMs));
    while($skip_random{$HMMs[$r]}){ $r = int(rand($#HMMs)) }
    $skip_random{$HMMs[$r]} = 1;
    print "# skip $HMMs[$r]\n";
  }
  
  printf("# total randomly skipped HMMs: %d\n\n",scalar(keys(%skip_random)));
}  
	
## if requested parse keggmap file including pathway names
my %pathways;

if($INP_keggmapfile)
{
	open(KEGGMAP,$INP_keggmapfile);
	while(<KEGGMAP>)
	{
		#pfam 	ko	pathway name
		#PF00890 K00394  1 Sulfite oxidation
    chomp;

    next if(/^PFAM/ || /^#/); # skip header and any commented lines

    my ($pf,$kegg,$path_number,$path_name) = split(/\t/,$_);
		push(@{$KEGGmap{$pf}},$kegg); 		
    
		$pathways{$path_number}{$pf} = 1;
    $pathways{$path_number}{'fullname'} = $path_name;
    $pathways{$path_number}{'totalHMMs'}++;
	}
	close(KEGGMAP);
}
	
## read input file with hmmsearch output in tab-separated format
my $maxmatches = 0;
my $pathwayOK;
if($INP_infile_bzipped)
{
	open(INFILE,"bzcat $INP_pfamsearchfile|") ||
		die "# $0 : cannot find $INP_pfamsearchfile, please check paths and re-run\n";
}
else
{
	open(INFILE,$INP_pfamsearchfile) ||
		die "# $0 : cannot find $INP_pfamsearchfile, please check paths and re-run\n";
}

while(<INFILE>)
{
	#5723145_1_98_+          -          2-Hacid_dh           PF00389.25   1.9e-08   36.7   0.1   1.9e-08   36.7   0.1   1.0   1   0   0   1   1   1   1 -
	#SRR000281.13791_1_100_+ -          2-Hacid_dh           PF00389.25   6.1e-06   28.6   0.1   6.1e-06   28.6   0.0   1.0   1   0   0   1   1   1   1 -
	chomp;
	next if(/^#/ || /^\s+$/);

	my @data = split(/\s+/,$_);
	$hmm = $data[3];  
	$hmm = (split(/\.\d+/,$hmm))[0];
	
	if($INP_pathways)
	{
		$pathwayOK = 0;
		foreach $pw (@user_pathways)
		{
			if($pathways{$pw}{$hmm}){ $pathwayOK = 1; last } 
		}
		next if($pathwayOK == 0);
	}
	
	if($HMMentropy{$hmm} && $HMMentropy{$hmm} >= $INP_minentropy)
	{
		$matchedHMMs{$hmm}++;
		if($matchedHMMs{$hmm} > $maxmatches){ $maxmatches = $matchedHMMs{$hmm} }
		#print "$hmm $HMMentropy{$hmm}\n";
	}
}
close(INFILE);


## produce final score based on observed matched HMMs and also
## compute pathway completness 
my ($score,$matches,$color,@mapscript,@pws,%previous,%completness) = (0,0);
if(keys(%pathways)){ @pws = sort {$a<=>$b} keys(%pathways) }

print "# Pfam\tentropy\t#matched_peptides\n";

# sort HMMs by #matches
foreach $hmm (sort {$matchedHMMs{$b} <=> $matchedHMMs{$a}} keys(%matchedHMMs))
{
  next if($skip_random{$hmm}); 

	$matches = $matchedHMMs{$hmm} || 0;
	$entropy = $HMMentropy{$hmm} || 0;
	
  # print raw data for individual Pfam HMM
  print "$hmm\t$entropy\t$matches\n";

	if($matches>0)
	{
    # compute completness of pathways where this domain participates
    foreach $pw (@pws)
    {
      next if(!$pathways{$pw}{$hmm});
      push(@{$completness{$pw}},$hmm);
    }

    # prepare KEGG mapping script
		$score += $entropy;
		$color = $COLORS[ int(($matches/$maxmatches)*$#COLORS) ];
		
		# prepare KEGG mappings for http://www.genome.jp/kegg-bin/show_pathway
		foreach $KEGGid (@{$KEGGmap{$hmm}})
		{
			next if($KEGGid eq '' || $previous{$KEGGid});

			push(@mapscript, "$KEGGid $color,black");
			$previous{$KEGGid}++; 
		}
	}
}

print "\n# Pfam entropy score: $score\n\n";

if($INP_keggmapfile)
{
  my $mean=0;

  print "# Pathway report\n";
  print "# path_number\tpath_name\ttotal_domains\tmatched\t%completeness\tmatched_Pfam_domains\n";
  foreach $pw (@pws)
  {
    my $pw_matched_Pfams = '';
    my $pw_matches = 0;
    my $pw_comp = 0;
    if($completness{$pw})
    { 
      $pw_matches = scalar(@{$completness{$pw}});
      $pw_matched_Pfams = join(',',@{$completness{$pw}});
      $pw_comp = 100*($pw_matches/$pathways{$pw}{'totalHMMs'});
    } 

    printf("%d\t%s\t%d\t%d\t%1.1f\t%s\n",
      $pw,$pathways{$pw}{'fullname'},
      $pathways{$pw}{'totalHMMs'},
      $pw_matches,
      $pw_comp,
      $pw_matched_Pfams);

    $mean+=$pw_comp;
  }
  
  printf("\n# mean pathway completeness: %1.1f\n",$mean/scalar(@pws));

  print "\n\n# Script to map these Pfam domains in KEGG->User Data Mapping.\n";
  print "# Colors are proportional to the number of Pfam matches and normalized with respect to the max.\n";
  print "# Note that a reference map must be selected first. For instance, Sulphur metabolism is:\n";
  print "# http://www.genome.jp/kegg-bin/show_pathway?map00920\n";
  print "# WARNING: note that several K numbers might map to the same protein.\n";
  print "# These cases are coloured in KEGG with the color of the last given K number.\n";
  print "# You might want to manually adjust the colors to correct this.\n";
  print join("\n",sort(@mapscript))."\n";
}
