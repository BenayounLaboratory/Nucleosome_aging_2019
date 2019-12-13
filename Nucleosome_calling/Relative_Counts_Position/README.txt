README

##2016-08-05

- Used GREAT to plot distance to TSS of 3vs29 1e-15 DANPOS calls
- plot gained vs. lost
- wc -l on bed files, then clean up in Excel

##2016-10-06

# get ALL nucleosome positions as background
cat /Volumes/MyBook_3/BD_aging_project/ChIP-seq/New_H3_lanes/Nucleosome_analysis/DANPOS/Liver_H3_positionning_aging_1e-30_2016-06-02/reference_positions.xls \
	 | perl -lane 'next if ($. == 1); print "$F[0]\t$F[1]\t$F[2]\tNucleosome_".($.-1)."\t"' > Liver_ALL_H3_nucleosome_positions.bed
	 
cat /Volumes/MyBook_3/BD_aging_project/ChIP-seq/New_H3_lanes/Nucleosome_analysis/DANPOS/NPCs_H3_positionning_aging_1e-30_2016-06-02/reference_positions.xls  \
	 | perl -lane 'next if ($. == 1); print "$F[0]\t$F[1]\t$F[2]\tNucleosome_".($.-1)."\t"' > NPCs_ALL_H3_nucleosome_positions.bed
	 
cat /Volumes/MyBook_3/BD_aging_project/ChIP-seq/New_H3_lanes/Nucleosome_analysis/DANPOS/Heart_H3_positionning_aging_1e-30_2016-06-02/reference_positions.xls  \
	 | perl -lane 'next if ($. == 1); print "$F[0]\t$F[1]\t$F[2]\tNucleosome_".($.-1)."\t"' > Heart_ALL_H3_nucleosome_positions.bed
	 
cat /Volumes/MyBook_3/BD_aging_project/ChIP-seq/New_H3_lanes/Nucleosome_analysis/DANPOS/OB_H3_positionning_aging_1e-30_2016-06-02/reference_positions.xls  \
	 | perl -lane 'next if ($. == 1); print "$F[0]\t$F[1]\t$F[2]\tNucleosome_".($.-1)."\t"' > OB_ALL_H3_nucleosome_positions.bed
	 
cat /Volumes/MyBook_3/BD_aging_project/ChIP-seq/New_H3_lanes/Nucleosome_analysis/DANPOS/Cerebellum_H3_positionning_aging_1e-30_2016-06-02/reference_positions.xls  \
	 | perl -lane 'next if ($. == 1); print "$F[0]\t$F[1]\t$F[2]\tNucleosome_".($.-1)."\t"' > Cerebellum_ALL_H3_nucleosome_positions.bed

====> did not work (only summit is there!)

will just use genomic background

##### 
# 2017-03-20
# rerun with corrected Cerebellum
