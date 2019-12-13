#! /usr/bin/perl

use warnings;
use strict;

# 2019-11-04
# parse HOMER Genome Ontology randomizations into one file

unless (@ARGV == 3) {
	die "\nparse_HOMER_GO.pl <outfile name> <HOMERGO Real> <HOMERGOfile Rand Folder>\n\n";
}

my $outname = shift @ARGV;
my $real = shift @ARGV;
my $randomsDIR = shift @ARGV;


#############################################################################
# open and parse real file
my %outresults = ();

open(REAL, $real) or die "Could not open $real: $!\n";

# get first line
my $headerline = <REAL>;
#Name	PeakFile/Annotation	#features[ref=11178]	Coverage(bp)[ref=1553457]	AvgFeatureSize[ref=138]	Overlap(#peaks)	Overlap(bp)	Expected Overlap(bp, gsize=2.00e+09)	Log Ratio Enrichment	Log P-value(+ underrepresented)	P-value

while (my $line = <REAL>) {

	# parse fields
	my @linedata = get_line_data($line);
	
	# $linedata[0] is the genomic element name
	# $linedata[1] is the homer annotation file
	# $linedata[2] is number of features
	# $linedata[3] is the genomic coverage of features
	# $linedata[4] is the AvgFeatureSize                             
	# $linedata[5] is the Overlap(#peaks)                      %#%#%# feature to extract                   
	# $linedata[6] is the Overlap(bp)                             
	# $linedata[7] is the Expected Overlap(bp, gsize=2.00e+09) #### HOMER values based on whole genome!!!!!
	# $linedata[8] is the Log Ratio Enrichment                 #### HOMER values based on whole genome!!!!!
	# $linedata[9] is Log P-value(+ underrepresented)          #### HOMER values based on whole genome!!!!!
	# $linedata[10] is P-value                                 #### HOMER values based on whole genome!!!!!
	push(@{$outresults{$linedata[0]}}, $linedata[5]);
}


close REAL;


#############################################################################
# read and parse files from random dir
opendir(my $dh, $randomsDIR) or die "Can't opendir $randomsDIR: $!";
my @randoms = readdir($dh);
closedir $dh;
##### print join("\t",@randoms)."\n";

my @randoms_Basic = grep(/_basic_genomeOntology.txt$/, @randoms);
#print join("\t",@randoms_Basic)."\n";
#print @randoms_Basic."\n";

# make output header
my $header = "Genomic_Element\tOveralp_Real_Number_Nucs\t";
$header .= join("\t",@randoms_Basic);
$header .= "\n";

# loop on files in random directory
foreach my $file (@randoms_Basic) {
	
	# add folder path to read file
	$file = $randomsDIR."/".$file;
	
	open(FILE,$file) or die "Could not open $file: $!\n";
    
    # skip first line
    my $headerline = <FILE>;
    
    # loop on lines
    while (my $line = <FILE>) {
		my @linedata = get_line_data($line);
		push(@{$outresults{$linedata[0]}}, $linedata[5] );
	}
    
    close FILE;
	
}


#############################################################################
#### output results
open(OUT,'>',$outname) or die "Could not open $outname: $!\n";

print OUT $header;

foreach my $elem (sort keys %outresults) {
	
	#print $elem."\n";
	
	my $outline = $elem."\t".join("\t",@{$outresults{$elem}})."\n";
	print OUT $outline;
}

close OUT;

exit;

###########################################################
# SUBROUTINES
###########################################################

###########################################################
# a subroutine that separates fields from a data line and
# returns them in an array

sub get_line_data {

    my $line = $_[0];
    
    chomp $line;  

    my @linedata = split(/\t/, $line);
        
    return @linedata;
}
