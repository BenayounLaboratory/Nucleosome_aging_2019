setwd('/Volumes/MyBook_3/BD_aging_project/ChIP-seq/All_tissues_analysis/Nucleosome_analyses/Relative_Counts//Consensus')

# 2016-11-18

# 2017-03-20
# rerun with corrected numbers (cerebellum, danpos reference position)


#  wc -l dposREF*DiNUP_DANPOS*bed
# 6557 dposREF_Cerebellum_DiNUP_DANPOS_GAINED.bed
# 6009 dposREF_Cerebellum_DiNUP_DANPOS_LOST.bed
# 3348 dposREF_Heart_DiNUP_DANPOS_GAINED.bed
# 4032 dposREF_Heart_DiNUP_DANPOS_LOST.bed
# 5759 dposREF_Liver_DiNUP_DANPOS_GAINED.bed
# 7777 dposREF_Liver_DiNUP_DANPOS_LOST.bed
# 2458 dposREF_NPCs_DiNUP_DANPOS_GAINED.bed
# 1412 dposREF_NPCs_DiNUP_DANPOS_LOST.bed
# 7950 dposREF_OlfactoryBulb_DiNUP_DANPOS_GAINED.bed
# 11178 dposREF_OlfactoryBulb_DiNUP_DANPOS_LOST.bed


my.nucs <- rbind(
                c(3348,5759,6557,7950,2458),
                c(4032,7777,6009,11178,1412)
                )
  

colnames(my.nucs) <- c("Heart","Liver","Cerebellum","OB","NPCs")
rownames(my.nucs) <- c("GainedNuc","LostNuc")

my.nucs.1 <- my.nucs
my.nucs.1[2,] <- -1* as.numeric(my.nucs.1[2,])

pdf("2017-03-20_barplot_changed_nucleosomes_consensus_DANPOS-DINUP.pdf", width=8, height=5)
par(oma=c(0.5,1.5,0.5,0.5))
barplot(as.matrix(my.nucs.1), beside=T, las=1,
        col = c("#CC3333","#333399"),horiz=T,
        xlab= "Significantly changed nucleosome positions (DANPOS and DINUP consensus)",
        cex.axis = 0.75, cex.names = 0.75, cex.lab= 0.75,
        xlim = c(-12000,12000))
abline(v=0, lty="dashed")
box()
dev.off()



