#! /usr/bin/perl

use warnings;
use strict;

my $matrix = "Nucleosome_data_matrix_v2.txt";

#my $out = "2016-08-25_joblist_for_GREAT_nucleosomes.txt";
my $out = "2016-11-17_joblist_for_GREAT_nucleosomes.txt";

# open filehandles
open (MATRIX, $matrix) or die "couldn't open $matrix: $!\n";
open (OUT, '>', $out) or die "couldn't create $out: $!\n";

my $header = <MATRIX>;
my $dropbox_folder = 'https%3A%2F%2Fdl.dropboxusercontent.com%2Fu%2F19880001%2F';

# create the parse
while (my $line = <MATRIX>) {
  #  http://bejerano.stanford.edu/great/public/cgi-bin/greatStart.php?
  #  outputType=batch&
  #  requestSpecies=hg18&
  #  requestName=Example+Data&
  #  requestSender=BereniceBenayoun&
  #  requestURL=http%3A%2F%2Fwww.clientA.com%2Fdata%2Fexample1.bed&
	#print $line;
	
	my ( $bedfile, $assembly, $prefix) = get_line_data($line);
	
	my $outline = "$prefix"."_GREAT.tsv\t".
				  "http://bejerano.stanford.edu/great/public/cgi-bin/greatStart.php?".
				  "outputType=batch&".
				  "requestSpecies=$assembly"."&".
				  "requestName=$prefix"."&".
				  "requestSender=BereniceBenayoun&".
				  "requestURL=$dropbox_folder".$bedfile
					;
					
	print OUT $outline."\n";
	

}


close MATRIX;
close OUT;

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
    
    my @linedata = split(/\s/, $line);
       
    return @linedata;
}
