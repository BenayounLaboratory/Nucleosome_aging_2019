#! /usr/bin/perl

use warnings;
use strict;
use File::Basename;

unless (scalar @ARGV == 3) {
	die "\nIncorrect command line arguments.\n".
		"Usage : get_random_nucleosomes.pl <repeats> <nuc number> <nucleosome background> \n";
}

my $repeats = shift @ARGV;
my $num     = shift @ARGV;
my $nucfile = shift @ARGV;

# read in nucleosome background
open (NUCS,$nucfile) or die "Couldn't open $nucfile: $!\n";
my @nucs = <NUCS>;
close NUCS;

my $infname = basename($nucfile,(".bed"));
my $oDir = $infname."_RANDOMS_".$repeats;
system("mkdir $oDir");

for (my $i =1; $i <= $repeats; ++$i) {
	
	my $oFname = "RANDOM_nucleosomes_".$i.".bed";
	
	# select num unique nucleosomes
	# modified from https://www.perlmonks.org/?node_id=244252
	my @rand_nucs = ();
	my %rand_nucs_hash = ();

	while (1) {
  		my $nuc = $nucs[ rand @nucs ];
  		
  		unless (defined $rand_nucs_hash{$nuc}) {
  			++$rand_nucs_hash{$nuc}; # make sure each original nucleosome is only used once
  			push (@rand_nucs, $nuc);
  		}
  		
  		last if (scalar @rand_nucs == $num);
	}
	
	
	open (OUT,'>',$oDir."/".$oFname) or die "Couldn't create $oFname: $!\n";
	print OUT @rand_nucs;
	close OUT;
}

exit;

#####################################
######      SUBROUTINES        ######
#####################################

