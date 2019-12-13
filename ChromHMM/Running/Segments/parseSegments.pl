#! /usr/bin/perl

use warnings;
use strict;

###############################################
# USAGE STATEMENT
unless (scalar @ARGV >= 2) {
	die "\nIncorrect command line arguments.\n".
		"Usage : parseSegments.pl <prefix> <HMM file> \n";
}
###############################################

my $prefix = shift @ARGV;
my $segmentfile = shift @ARGV;

my %segment_hash = ();


open (SEGMENTS,$segmentfile) or die "Couldn't open $segmentfile: $!\n";

while (my $line = <SEGMENTS>) {
	my @linedata = get_line_data ($line);
	
	push(@{$segment_hash{$linedata[3]}},$line);
	
}

close SEGMENTS;



foreach my $state (sort keys %segment_hash) {
	
	my $outfile = $prefix."_ChromHMM_".$state.".bed";
	
	open (BED,'>',$outfile) or die "Couldn't open $outfile: $!\n";
	
	print BED @{$segment_hash{$state}};
	close BED;

}


exit;

#####################################
######      SUBROUTINES        ######
#####################################

###########################################################
# a subroutine that separates fields from a data line and
# returns them in an array

sub get_line_data {

    my $line = $_[0];
    
    chomp $line;
    
    my @linedata = split(/\t/, $line);
       
    return @linedata;
}
