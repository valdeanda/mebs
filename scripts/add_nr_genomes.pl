#!/usr/bin/perl

use strict;

die "usage: $0 <list of nr genome uids> <folder with NCBI genome<>\n" if(!$ARGV[1]);

my ($NRLIST,$NCBI_RAW_FOLDERS) = @ARGV;

my (%nr,$uid);
open(NRLIST,$NRLIST);
while(<NRLIST>)
{
    #uid156845,
    next if(/^#/);
    chomp;
    $uid = (split(/,/,$_))[0]; #print "|$uid|\n";
    $nr{$uid}++;
}
close(NRLIST);

#printf("# total nr genomes in list = %d\n",scalar(keys(%nr)));

opendir(RAW,$NCBI_RAW_FOLDERS);
my @genomedirs = grep {/uid/} readdir(RAW);
closedir(RAW);

foreach my $dir (@genomedirs)
{
    if($dir =~ /(uid\d+)/)
    {
        $uid = $1; #print "$uid|\n";
        next if(!$nr{$uid});
        #print "$dir\n";

        opendir(UIDDIR,"$NCBI_RAW_FOLDERS/$dir");
        my @faafiles = grep{/faa/} readdir(UIDDIR);
        closedir(UIDDIR);

        foreach my $faafile (@faafiles)
        {
            open(FAA,"$NCBI_RAW_FOLDERS/$dir/$faafile");
            while(<FAA>)
            {
                print;
            }
            close(FAAFILE);
        }
    }   
}

