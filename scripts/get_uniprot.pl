#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use FindBin '$Bin';
#Modify from https://www.uniprot.org/help/api_idmapping#id_mapping_perl_example
use LWP::UserAgent;

my ($INP_help, $INP_file, $INP_opt) = (0,'','');
my $INP_default = "BIOCYC_ID";
#UNIPATHWAY_ID
#BIOCYC_ID


GetOptions
(
  'help|h|?'    => \$INP_help,
  'input|in=s'  => \$INP_file,
  'option|s:o'  => \$INP_opt,
);



if($INP_help || $INP_file eq '')
{
  die<<EODOC;

  Program to Mapping database identifiers.

  usage: $0 [options] 

   -help    Brief help message
   
   -input   File containing Uniprot identifiers  (required)

   -option  Database to retrieve identifier      (required, default $INP_default)
            see https://www.uniprot.org/help/api_idmapping#id_mapping_perl_example
             
EODOC
}


#Required input files  

if (!$INP_file || !-s $INP_file)
{
  die  "# Error: need valid -input $INP_file \n";
}

print "# $0 call:\n# -input $INP_file -option $INP_opt\n\n";




my (@identifiers);
#Read input file  

open (INFILE,$INP_file)  || die "# $0,  ERROR: cannot read $INP_file\n";
while (my $line =<INFILE>) 
{
  chomp $line;
  my @file = split(/\t/,$line);
  push (@identifiers,$file[0]);
}
close (INFILE);

#print "# Your input identifiers are:\n".join ("\n",@identifiers)."\n\n"; 

if (!$INP_opt)
{
my $base = 'http://www.uniprot.org';
my $tool = 'uploadlists';

my $params = {
  from => 'ACC+ID',
  to => $INP_default,
  format => 'tab',
  query  => join(' ',@identifiers),
};

my $contact = ''; # Please set your email address here to help us debug in case of problems.
my $agent = LWP::UserAgent->new(agent => "libwww-perl $contact");
push @{$agent->requests_redirectable}, 'POST';

my $response = $agent->post("$base/$tool/", $params);

while (my $wait = $response->header('Retry-After')) {
  print STDERR "Waiting ($wait)...\n";
  sleep $wait;
  $response = $agent->get($response->base);
}
#print "Identified IDs\n";
$response->is_success ?
  print $response->content :
  die 'Failed, got ' . $response->status_line .
    ' for ' . $response->request->uri . "\n";
}

if ($INP_opt) 

{
 my $base = 'http://www.uniprot.org';
 my $tool = 'uploadlists';
 my $params = {
  from => 'ACC+ID',
  to => $INP_opt,
  format => 'tab',
  query  => join(' ',@identifiers),
};

  my $contact = ''; # Please set your email address here to help us debug in case of problems.
  my $agent = LWP::UserAgent->new(agent => "libwww-perl $contact");
  push @{$agent->requests_redirectable}, 'POST';
  my $response = $agent->post("$base/$tool/", $params);
  while (my $wait = $response->header('Retry-After'))
  {
    print STDERR "Waiting ($wait)...\n";
    sleep $wait;
    $response = $agent->get($response->base);
   }
   #print "Identified IDs\n";
   $response->is_success ?
   print $response->content :
   die 'Failed, got ' . $response->status_line .
    ' for ' . $response->request->uri . "\n";
 }









