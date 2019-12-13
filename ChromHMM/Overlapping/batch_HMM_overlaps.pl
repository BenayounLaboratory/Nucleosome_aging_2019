#! /usr/bin/perl

use warnings;
use strict;
use IO::CaptureOutput qw(capture_exec);

unless (scalar @ARGV >= 2) {
	die "\nIncorrect command line arguments.\n".
		"Usage : batch_HMM_overlaps.pl <CellType> <NucFileUP> <NucFileDWN> \n";
}

my $CellType = shift @ARGV;
my $gained = shift @ARGV;
my $lostnuc = shift @ARGV;

my $segDir = '/Volumes/MyBook_3/BD_aging_project/Public_datasets/LICR_Datasets/ChromHMM/CHROM_HMM_RUN_v3/Segments';

# get segments
opendir(my $dh, $segDir) || die "Can't opendir $segDir: $!";
my @Segments = grep { /^$CellType/ && !/HMM_segments\.bed/ } readdir($dh);
closedir $dh;

#my $segnums = scalar @Segments;
#print  $segnums."\n\n";
#print join("\n",@Segments)."\n\n";

#create report file
my $report = $CellType."_nucleosomes_HMM_report.txt";

open (FILE,'>',$report) or die "Couldn't open $report: $!\n";
print FILE "$CellType\n";


#### 1. gained nucs
# number of nucleosomes
my $numTOT = capture_exec("wc -l $gained");
$numTOT =~ s/\s+/\t/g; $numTOT =~ s/^\t//g;
my @numTOT = split(/\s/,$numTOT);

#print join("::",@numTOT)."\n\n";

print FILE "$numTOT[0] Gained Nucleosomes\n";

for (my $i = 0; $i < scalar @Segments; ++$i) {
	my $file = $Segments[$i];

	$file =~ m/_ChromHMM_(E\d+)\.bed/;
	my $state = $1;
	my $outname = $state."_".$CellType."_gained_nucleosomes.bed";
	my $cmd = "intersectBed -f 0.50 -a $gained -b $segDir/$file > $outname";
	system($cmd);
	
	# now get number of nucleosomes
	my $num = capture_exec("wc -l $outname");
	$num =~ s/\s+/\t/g; $num =~ s/^\t//g;
	my @num = split(/\s/,$num);

	print FILE "$state\t$num[0]\n";

}



###### 2. lost nucs
# number of nucleosomes
$numTOT = capture_exec("wc -l $lostnuc");
$numTOT =~ s/\s+/\t/g; $numTOT =~ s/^\t//g;
@numTOT = split(/\s/,$numTOT);

#print join("\n",@Segments)."\n\n";

print FILE "\n\n$numTOT[0] Lost Nucleosomes\n";

for (my $i = 0; $i < scalar @Segments; ++$i) {
	my $file = $Segments[$i];

	$file =~ m/_ChromHMM_(E\d+)\.bed/;
	my $state = $1;
	my $outname = $state."_".$CellType."_lost_nucleosomes.bed";
	my $cmd = "intersectBed -f 0.50 -a $lostnuc -b $segDir/$file > $outname";
	system($cmd);
	
	# now get number of nucleosomes
	my $num = capture_exec("wc -l $outname");
	$num =~ s/\s+/\t/g; $num =~ s/^\t//g;
	my @num = split(/\s/,$num);

	print FILE "$state\t$num[0]\n";

}


close FILE;

exit;

#####################################
######      SUBROUTINES        ######
#####################################

