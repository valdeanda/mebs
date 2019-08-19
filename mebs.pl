#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use FindBin '$Bin';

# General script to score biogeochem cycles in both genomic and metagenomic data.
# B Contreras-Moreira, V de Anda 2018
# bcontreras@eead.csic.es , valdeanda@ciencias.unam.mx 

my $VERSION = 'v1.0';
my $DEBUG = 0;

my $HMMSEARCHEXE = 'hmmsearch'; # please edit if not in path

my $CONFIGDIR   = $Bin.'/config/';
my $CONFIGFILE  = $CONFIGDIR . 'config.txt';
my $VALIDEXT    = '.faa';
my $VALIDENT    = '.tab'; # valid extension for pre-computed entropy files
my $VALIDHMMEXT = '.hmm'; # valid extension for HMM files
my $HMMOUTEXT   = '.hmmsearch.tab'; # default extension for hmmsearch tbl output
my $FDR         = 0.01;
my @validFDR    = qw(0.1 0.01 0.001 0.0001);
my @validMSL    = qw(30 60 100 150 200 250 300);
my ($INP_help,$INP_folder,$INP_cycles,$INP_type,$INP_FDR, $INP_comp) = (0,'',0,'',$FDR,0);

GetOptions
(
  'help|h|?'    => \$INP_help,
  'input|in=s'  => \$INP_folder,
  'type|t=s'    => \$INP_type,
  'cycles|c'    => \$INP_cycles,
  'fdr|r=f'     => \$INP_FDR,
  'comp|mc'    => \$INP_comp,
);


if (-t STDIN && ($INP_help || $INP_folder eq '' || $INP_type eq '') && !$INP_cycles)
{
  die<<EODOC;

  Program to compute MEBS for a set of genomic/metagenomic FASTA files in input folder.
  Version: $VERSION

  usage: $0 [options] 

   -help    Brief help message
   
   -input   Folder containing FASTA peptide files ($VALIDEXT)             (required)

   -type    Nature of input sequences, either 'genomic' or 'metagenomic'  (required)

   -fdr     Score cycles with False Discovery Rate @validFDR  (optional, default=$FDR)

   -cycles  Show currently supported biogeochemical cycles/pathways
   
   -comp    Compute the metabolic completeness of default cycles.         (optional) 
            Required option  for mebs_output.py                                 

   -custom  Compute the metabolic completeness of user input pathways     (optional) 

EODOC
}

## 1) Checking binaries
my $sample_output = `$HMMSEARCHEXE -h 2>&1 `;
if(!$HMMSEARCHEXE || $sample_output !~ /Usage/)
{
  die "#ERROR:  hmmsearch not found, please install or set \$HMMSEARCHEXE correctly\n";
}

## 2) Checking parameters
my (@valid_infiles, @cycles, @config, @paths, @completeness);
my (@MSL, %FDRcutoff, %col2fdr, %pathways);
my ($c,$f,$path,$cycle,$msl,$score,$comp,$pw);
my ($hmmfile,$hmmsearchfile,$entropyfile,$scorefile,$infile,$pfam2keggfile);


# Read config file
open(CONFIG,$CONFIGFILE) || die "# ERROR: cannot read $CONFIGFILE\n";
while(my $line = <CONFIG>)
{
  #Cycle   Path    Comple  Input Genes     Input Genomes   Domains AUC     Score(FD..
  #sulfur  cycles/sulfur/  cycles/sulfur/pfam2kegg.tab     152     161     112    ..
  #oxygen  cycles/oxygen/    50  53  55  ...

  next if($line =~ m/^\s+/);

  @config = split(/\t/,$line);
  if($config[0] =~ /^Cycle/)
  {
    # check which columns in config match which FDR-based cutoffs
    foreach $c (1 .. $#config)
    {
      if($config[$c] =~ /Score\(FDR(0\.\d+)/)
      {
        $col2fdr{$c} = $1;
      }
    }
  }
  else
  {
    if($DEBUG == 1){ print "$config[0],$config[1],$config[2]\n" }
    push(@cycles, $config[0]);
    push(@paths, $config[1]);
    push(@completeness, $config[2]); # $config[2] might be '' if not curated 


    # save score FDR cutoffs
    foreach $c (keys(%col2fdr))
    {
      $FDRcutoff{$config[0]}{$col2fdr{$c}} = $config[$c];
    }
  }  
}
close(CONFIG); 


if ($INP_cycles)
{
  print "# Available cycles:\n". join("\n",@cycles)."\n\n";
  print "# Available files to compute completeness:\n". join("\n",@completeness)."\n\n";
  exit(0);
}
else
{
  warn "# $0 -input $INP_folder -type $INP_type -fdr $INP_FDR -comp $INP_comp\n\n";
}
 
# check required sequence type
if(!$INP_type || ($INP_type ne 'genomic' && $INP_type ne 'metagenomic'))
{
  die "# ERROR : type of input must be indicated; valid options are [genomic|metagenomic]\n";
}

# check required sequence folder
if(!$INP_folder)
{
  die "# ERROR : need valid -input folder\n";
}
else
{
  opendir(DIR,$INP_folder) || die "# ERROR: cannot list $INP_folder\n";
  @valid_infiles = grep{/$VALIDEXT$/} readdir(DIR);
  closedir(DIR);
  if(scalar(@valid_infiles) == 0)
  {
    die "# ERROR: cannot find files with extension $VALIDEXT in folder $INP_folder\n";
  }
  elsif($INP_type eq 'metagenomic')
  {
    # compute Mean Size Length for this metagenomic sequence set
    warn "# Computing Mean Size Length (MSL) ...\n";

    my ($c,$nseq,$mean,$len,$cutoff,@MSLcutoffs);
    for(my $bin=0;$bin<scalar(@validMSL)-1;$bin++)
    {
      $cutoff = $validMSL[$bin] + (($validMSL[$bin+1]-$validMSL[$bin])/2);
      push(@MSLcutoffs,$cutoff);#print "$validMSL[$bin] $cutoff\n";
    }
    push(@MSLcutoffs,$validMSL[$#validMSL]); # add last MSL

    foreach my $infile (@valid_infiles)
    {
      ($nseq,$mean,$len) = (0,0,0);
      open(FAAFILE,"<","$INP_folder/$infile") || 
        die "# ERROR: cannot read files $INP_folder/$infile\n";
      while(<FAAFILE>)
      {
        if(/^>/)
        {
          $nseq++;
          $mean += $len;
          $len=0;
        }
        else
        {
          chomp;
          $len += length($_);
        }
      }
      close(FAAFILE);

      $mean = sprintf("%1.0f",$mean/$nseq);

      # find out which pre-defined MSL bin matches the estimated MSL for this sequence set
      foreach $c (0 .. $#MSLcutoffs)
      {
        $cutoff = $MSLcutoffs[$c];
        if($mean <= $cutoff)
        {
          push(@MSL,$validMSL[$c]);
          warn "# $infile MSL=$mean MSLbin=$validMSL[$c]\n";

          last;
        }
      }
    }
  } print "\n";
}

# check optional FDR 
if($INP_FDR)   
{
  if(!grep (/^$INP_FDR$/, @validFDR))
  {
    die "# ERROR: FDR value is not valid; please choose from ".join(', ',@validFDR)."\n";
  }
}

###Option to use TGRFAM or Pfam 

##if custom is selected download current pfam database  
#option would you like to dowload current Pfam database?
##Warning heavy file 1.4 G
#Download PFAM => ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/
#gzip Pfam-A.hmm.gz
#mv Pfam-A.hmm.gz my_Pfam.custom.hmm && mv my_Pfam.custom.hmm cycles/custom
#
#
#
## 3) scan input sequences with selected Pfam HMMs for each input file & cycle

# print header
my $pathways_header = '';
foreach $c (0 .. $#cycles)
{
  print "\t$cycles[$c]";

  # print completeness header if required
  $comp = $completeness[$c] || "";  
  if($INP_comp && $comp ne "" && -s $comp)
  {
    open(COMPFILE,"<",$comp) || warn "# ERROR: cannot read $comp\n";
    while(<COMPFILE>)
    {
      #PFAM  KO  PATHWAY   PATHWAY NAME 
      #PF00890 K00394  1 Sulfite oxidation 
      if(/^PF\d+\t.*?\t(\d+)\t/)
      {
        $pathways{$cycles[$c]}{$1} = 1; 
      }
    }
    close(COMPFILE);
     
    $pathways_header .="\t<$cycles[$c] comp>";
    foreach $pw (sort {$a<=>$b} keys(%{$pathways{$cycles[$c]}}))
    {
      $pathways_header .= "\t$cycles[$c]_$pw";
    }
  }
} print "$pathways_header\n"; 

foreach $f (0 .. $#valid_infiles)
{
  $infile = $valid_infiles[$f];

  print "$infile"; # rowname

  # compute & print scores per cycle
  foreach $c (0 .. $#cycles)
  {
    $path = $paths[$c];
    $cycle = $cycles[$c];
    $comp = $completeness[$c];
    $score = 'NA';

    $hmmsearchfile = $INP_folder . '/' . $infile . '.' . $cycle . $HMMOUTEXT;
    $scorefile = $INP_folder . '/' . $infile . '.' . $cycle . '.score';
    $hmmfile = $path . 'my_Pfam.'. $cycle . $VALIDHMMEXT;
    $entropyfile = $path . 'entropies' . $VALIDENT;

    if(!-s $hmmsearchfile)
    {
      system("$HMMSEARCHEXE --cut_ga -o /dev/null --tblout $hmmsearchfile $hmmfile $INP_folder/$infile");
    }

    if(-s $hmmsearchfile)
    {
      if($INP_type eq 'metagenomic')
      {
        if($INP_comp && $comp ne "" && -s $comp)
        { 
          system("$Bin/scripts/pfam_score.pl -input $hmmsearchfile -entropyfile $entropyfile -size $MSL[$f] -keggmap $comp > $scorefile");
        }
        else
        {
          system("$Bin/scripts/pfam_score.pl -input $hmmsearchfile -entropyfile $entropyfile -size $MSL[$f] > $scorefile");
        }
      }
      else
      {
        if ($INP_comp && $comp ne "" && -s $comp)
        {
          system("$Bin/scripts/pfam_score.pl -input $hmmsearchfile -entropyfile $entropyfile -size real -keggmap $comp > $scorefile");
        }
        else
        {
          system("$Bin/scripts/pfam_score.pl -input $hmmsearchfile -entropyfile $entropyfile -size real > $scorefile");
        }  
      }
      
      if(-s $scorefile)
      {
        $score = -1;
        open(SCORE,"<",$scorefile) || warn "# ERROR: cannot read $scorefile\n";
        while(<SCORE>)
        {
          if(/Pfam entropy score: (\S+)/){ $score = sprintf("%1.3f",$1) } 
        }
        close(SCORE);  
      }
      else { warn "# ERROR: failed to generate $scorefile\n"  }
    }
    else { warn "# ERROR: failed to generate $hmmsearchfile\n" }

    # compare score to FDR-based cutoff
    if($score ne 'NA' && $score >= $FDRcutoff{$cycle}{$INP_FDR})
    {
      $score .= '*';
    }

    print "\t$score";
  }

  # print completeness summary per cycle
  if ($INP_comp)
  {
    foreach $c (0 .. $#cycles)
    {
      $cycle = $cycles[$c];

      # parse score file and check whether completeness report is there
      my ($compOK,%comp_scores) = (0);
      $scorefile = $INP_folder . '/' . $infile . '.' . $cycle . '.score';
      open(COMPL,"<",$scorefile);
      while(<COMPL>)
      {   
        # path_number path_name total_domains matched %completeness matched_Pfam_domains
        #1 Sulfite oxidation   9 3 33.3  PF00890,PF12838,PF01087
        # mean pathway completeness: 37.9
        if(/^# path_number/){ $compOK = 1 }
        elsif($compOK && /^(\d+)\t.*?\t\d+\t\d+\t(\S+)/)
        {
          $comp_scores{$1} = $2;
        }
        elsif(/^# mean pathway completeness: (\S+)/)
        {
          print "\t$1"; # print mean 
        }
      }
      close(COMPL);
    
      # print completeness for all sorted pathways
      foreach $pw (sort {$a<=>$b} keys(%{$pathways{$cycle}}))
      {
        print "\t$comp_scores{$pw}";
      }  
    }
  }    


  print "\n";
}























