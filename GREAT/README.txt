# 2016-08-25

Use GREAT programmatic interface to get data (greatBatchQuery.py)
# http://bejerano.stanford.edu/help/display/GREAT/Programming+Interface
# http://davetang.org/muse/2014/05/16/genomic-regions-enrichment-of-annotations-tool/

#Create a joblist of the format: outputFile queryURL where each line is a new job
#Then run ./greatBatchQuery.py joblist


copy files into the Dropbox public folder to get the URL for batch submission
https://dl.dropboxusercontent.com/u/19880001/


./greatBatchQuery.py 2016-08-25_joblist_for_GREAT.txt

get_broad_from_bed.pl 0 Heart_background_H3K4me3_breadth_ForGREAT.bed Lens_background_H3K4me3_breadth_ForGREAT.bed Liver_background_H3K4me3_breadth_ForGREAT.bed MuSCs_background_H3K4me3_breadth_ForGREAT.bed NPCs_background_H3K4me3_breadth_ForGREAT.bed OB_background_H3K4me3_breadth_ForGREAT.bed PFC_neurons_background_H3K4me3_breadth_ForGREAT.bed DF_background_H3K4me3_breadth_ForGREAT.bed HSCs_background_H3K4me3_breadth_ForGREAT.bed Cerebellum_background_H3K4me3_breadth_ForGREAT.bed


remove manually headers in background
need to concatenate foreground and background for GREAT to run

#bedEndRepair.pl /Users/benayoun/Softwares/BedTools-2.16.1/genomes/mouse.mm9.genome Olfactory_Bulb_3m1_H3K4me3_ADJ_broad_peaks.bed


# 2016-11-17
# rerun with consensus nucleosomes
./greatBatchQuery.py ./ForGREAT_nucleosomes/2016-11-17_joblist_for_GREAT_nucleosomes.txt
