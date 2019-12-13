setwd('/Volumes/MyBook_3/BD_aging_project/Public_datasets/LICR_Datasets/ChromHMM/CHROM_HMM_RUN_v3/Segments')


my.files <- list.files(pattern="bed")

my.results <- rep(0,length(my.files))

for ( i in 1:length(my.files)) {
  
  my.data <- read.table(my.files[i],sep="\t",header=F)
  
  my.results[i] <- sum(as.numeric(my.data$V3-my.data$V2))
  
}

write.table(cbind(my.files,my.results), file = "2016-12-22_Genome_occupancy_basepairs_segments.txt", quote=F,sep="\t",col.names=F, row.names=F)



# http://genomewiki.ucsc.edu/index.php/Genome_size_statistics
# mm9 : 2,725,765,481bp

my.results/2725765481
write.table(cbind(my.files,round(my.results/2725765481, digits=5)), file = "2016-12-22_Genome_occupancy_percentMM9_segments.txt", quote=F,sep="\t",col.names=F, row.names=F)

pdf('barplot_genome_occupancy.pdf', width=20, height=7.5)
par(oma=c(1,15,1,1))
barplot(100*my.results/2725765481, names.arg=my.files, horiz=T, las=2, xlab = "% mm9 genome")
dev.off()

pdf('barplot_genome_occupancy_log_scale.pdf', width=20, height=7.5)
par(oma=c(1,15,1,1))
barplot(100*my.results/2725765481, names.arg=my.files, horiz=T, las=2, xlab = "% mm9 genome", log = 'x', col = "tomato")
box()
dev.off()
