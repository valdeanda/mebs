#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $DEFAULTFRAGSIZE = 100;
my $DEFAULTCOVER    = 1;

my ($INP_infaafile,$INP_outfaafile) = ('','');
my ($INP_fragment_size,$INP_cover,$INP_help) = ($DEFAULTFRAGSIZE,$DEFAULTCOVER,0);

GetOptions
(
    'help|h|?'        => \$INP_help,
    'inFASTA|in=s'    => \$INP_infaafile,
	 'outFASTA|out=s'  => \$INP_outfaafile,
    'size|s:i'        => \$INP_fragment_size,
    'cover|c:i'       => \$INP_cover,
);

if (-t STDIN && ($INP_help || $INP_infaafile eq '' || $INP_outfaafile eq ''))
{
die<<EODOC;

Program to produce random fragments of proteins in input file with size and coverage set by user.

usage: $0 [options] 

 -help             brief help message
 
 -inFASTA          input file with protein sequences in FASTA format

 -outFASTA         output file with protein fragments in FASTA format
 
 -size             desired size for produced random fragments    (integer, default $INP_fragment_size)
 
 -cover            desired protein coverage of produced fragment (integer, default $INP_cover) 
						 
EODOC
}
   
if(!-s $INP_infaafile)
{
    die "# ERROR : cannot locate input file -fna $INP_infaafile\n";
}  

if($INP_outfaafile eq '')
{
	die "# ERROR : please provide a name for output file\n";
}

if($INP_cover < 1)
{
    die "# ERROR : invalid value for coverage ($INP_cover)\n";
}

if($INP_fragment_size < 1)
{
    die "# ERROR : invalid value for fragment size ($INP_fragment_size)\n";
}

print "# $0 call:\n# -inFASTA $INP_infaafile -outFASTA $INP_outfaafile -size $INP_fragment_size -cover $INP_cover\n\n";


my ($header,$seq,$n_of_frags,$coord,$length,$frag_length,%in_sequences);

# read input file
open(INFAA,$INP_infaafile);
while(<INFAA>)
{
	chomp;
	next if(/^#/ || /^\s+$/);
	
	if(/^>(.*)/){ $header = $1 }
	else{ $in_sequences{$header} .= $_; }
}
close(INFAA);

# loop through input sequences, produce fragments and save them in outfile

open(OUTFAA,">$INP_outfaafile");
foreach $header (keys(%in_sequences))
{
	$length = length($in_sequences{$header});
	$frag_length = $INP_fragment_size;
	$n_of_frags = 0;
	while($n_of_frags < int($INP_cover))
	{
		$coord = int(rand($length - $frag_length - 1));
		if($coord < 0){ $coord = 0; $frag_length = $length; }
		
		printf OUTFAA (">%s <fragment:%d,%d>\n",$header,$coord+1,$frag_length);
		print OUTFAA substr($in_sequences{$header},$coord,$frag_length)."\n";
		$n_of_frags++;
	}
}
close(OUTFAA);
