use strict;

my $IDLIST = 'id_interpro.txt';

my %DBs = ( 
  'Pfam'=>'Pfam-A.hmm',
	#'SUPERFAMILY'=>'interproscan-5.4-47.0/data/superfamily/1.75/hmmlib_1.75',
	#'TIGRFAM'=>'interproscan-5.4-47.0/data/tigrfam/13.0/TIGRFAMs_13.0_HMM.LIB'
);

my (%IDs,$db,$name,$hmmtext,$n_of_hmms);
open(IDLIST,$IDLIST);
while(<IDLIST>)
{
	#DATABASE	ID_DB
	#TIGRFAM	TIGR03567
	#SUPERFAMILY	SSF52218
	chomp $_;
	($db,$name) = split(/\t/,$_);
	if($db eq 'SUPERFAMILY')
	{
		$name =~ s/SSF//;
		#$name = sprintf("%07d",$name);
	}
	#print "|$db|$name|\n";
	$IDs{$db}{$name} = 1;	
}
close(IDLIST);  

foreach $db (keys(%DBs))
{
	print "# $db\n";
	
	open(MYDB,">my_$db.hmm");
	
	$n_of_hmms = 0;
	$hmmtext = $name = '';
	open(DB,$DBs{$db});
	while(<DB>)
	{
		if(/^\/\//)
		{
			$hmmtext .= $_;
			if($IDs{$db}{$name})
			{
				#print ">> $name\n";
				print MYDB $hmmtext;
				$n_of_hmms++ if($IDs{$db}{$name} == 1);
				$IDs{$db}{$name}++;
			}
			$hmmtext = '';
		}
		else
		{
			$hmmtext .= $_;	
			if($db =~/TIGRFAM/ && /^NAME\s+(\S+)/)
			{
				$name = $1;
				#print ">$name|\n";
			}
			elsif($db =~ /SUPERFAMILY/ && /^ACC\s+(\S+)/)
                        {
                                $name = $1;
				#print ">$name|\n";
                        }
			elsif($db =~/Pfam/  && /^ACC\s+(PF\d+)/)
                        {
                                $name = $1;
                                #print ">$name|\n";
                        }
		}
	}
	close(DB);	

	close(MYBD);

	print "# hmms = $n_of_hmms\n";

	foreach $name (keys(%{$IDs{$db}}))
	{
		print "falta $name\n" if ($IDs{$db}{$name}<2);
	}
}

	
