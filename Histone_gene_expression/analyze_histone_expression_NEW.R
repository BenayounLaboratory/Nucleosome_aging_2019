setwd('/Volumes/BB_Backup_3/BD_aging_project/RNAseq/All_tissues_analysis/Heatmaps/Histone_expression/')

# 2019-07-23
# look at histone gene expression with age


library('pheatmap')
library('beeswarm')

# read RNA log2 count matrix, DESeq2 normalized
my.rna.all <- read.csv('/Volumes/BB_Backup_3/BD_aging_project/RNAseq/All_tissues_analysis/DEseq2_runs/Global_and_nested/with_NPCs/2016-02-05_ALL_global_variance estimate_DESeq2_LINEAR_model_with_age _log2_counts_matrix.txt',
                       header=T,sep="\t")
my.rna.all$GeneName <- rownames(my.rna.all)

# set indexes for each
my.liv.ix <- 1:9
my.heart.ix <- 10:18
my.cereb.ix <- 19:27
my.ob.ix <- 28:35
my.npc.ix <- 36:41

my.histones <- read.csv('2019-07-23_Histone_Genes_Types.txt',header=T,sep="\t")
hist.idx <- which(rownames(my.rna.all) %in% my.histones$Associated.Gene.Name)

my.histone.exp <- merge(my.rna.all[hist.idx,], my.histones, by.x = "GeneName", by.y = "Associated.Gene.Name")
rownames(my.histone.exp) <- my.histone.exp$GeneName

########################################################################################################
# my.table   <- my.rna.all[hist.idx,my.liv.ix]

get_hist_norm <- function(my.table) {
  
  my.table <- my.table[apply(my.table, 1, sum)>0,]
  
  my.age <- rep("3m",ncol(my.table))
  my.age[grep("12m",colnames(my.table))] <- "12m"
  my.age[grep("29m",colnames(my.table))] <- "29m"
  
  # normalize to young level
  my.table.t <- data.frame(t(my.table))
  my.table.t <- 2^my.table.t
  my.mean.young <- apply(my.table.t[my.age %in% "3m",],2,mean)
  
  my.table.t.norm <- my.table.t
  for (t in 1:ncol(my.table.t)) {
    my.table.t.norm[,t] <- my.table.t[,t]/my.mean.young[t]
  }
  
  my.table.t.norm <- log2(my.table.t.norm)
  
  return(my.table.t.norm)
}



get_hist_plot <- function (my.liv.hist) {
  par(mfrow = c(1,ncol(my.liv.hist)))
  par(cex = 0.6)
  par(mar = c(0, 0, 0, 0), oma = c(4, 4, 10, 1))
  par(tcl = -0.25)
  par(mgp = c(2, 0.6, 0))
  
  for (i in 1:ncol(my.liv.hist)) {
    my.gene <- colnames(my.liv.hist)[i]
    
    my.age <- rep("3m",nrow(my.liv.hist))
    my.age[grep("12m",rownames(my.liv.hist))] <- "12m"
    my.age[grep("29m",rownames(my.liv.hist))] <- "29m"
    
    my.data.hist <- data.frame("gene" = my.liv.hist[,i],
                               "age" = my.age)
    colnames(my.data.hist)[1] <- "gene"
    
    
    boxplot(gene ~ age, data = my.data.hist,col=c("dodgerblue","purple","coral"),cex=2,axes = FALSE, ylim = c(-8,8))
    beeswarm(gene ~ age, data = my.data.hist,cex=1,col="black", pch = 16, add = TRUE)
    
    mtext(my.gene, side = 1, line = -2, adj = 0.1, cex = 1, col = "black", las =2)
    abline(h = 0, col = "red", lty = "dashed")
    
    if (i == 1)
      axis(2, col = "black", col.axis = "black", at = seq(-8,8,2),las=2)
    box(col = "black")
  }
  
  mtext("Relative expression by RNA-seq normalized to 3m", side = 2, outer = TRUE, line = 2.2,col = "black")
  
  par(mfrow = c(1,1))
  
}


########################################################################################################

# # select relevant rows
# my.liv.hist   <- get_hist_norm(my.histone.exp[hist.idx,my.liv.ix])
# my.heart.hist <- get_hist_norm(my.histone.exp[hist.idx,my.heart.ix])
# my.cereb.hist <- get_hist_norm(my.histone.exp[hist.idx,my.cereb.ix])
# my.ob.hist    <- get_hist_norm(my.histone.exp[hist.idx,my.ob.ix])
# my.npc.hist   <- get_hist_norm(my.histone.exp[hist.idx,my.npc.ix])
# 

# select relevant rows
my.liv.hist.H3   <- get_hist_norm(my.histone.exp[my.histone.exp$Histone.type %in% "H3",1+my.liv.ix]   )
my.heart.hist.H3 <- get_hist_norm(my.histone.exp[my.histone.exp$Histone.type %in% "H3",1+my.heart.ix] )
my.cereb.hist.H3 <- get_hist_norm(my.histone.exp[my.histone.exp$Histone.type %in% "H3",1+my.cereb.ix] )
my.ob.hist.H3    <- get_hist_norm(my.histone.exp[my.histone.exp$Histone.type %in% "H3",1+my.ob.ix]    )
my.npc.hist.H3   <- get_hist_norm(my.histone.exp[my.histone.exp$Histone.type %in% "H3",1+my.npc.ix]   )


pdf("Liver_H3_expression.pdf")
get_hist_plot(my.liv.hist.H3)
dev.off()

pdf("Heart_H3_expression.pdf")
get_hist_plot(my.heart.hist.H3)
dev.off()

pdf("Cereb_H3_expression.pdf")
get_hist_plot(my.cereb.hist.H3)
dev.off()

pdf("OB_H3_expression.pdf")
get_hist_plot(my.ob.hist.H3)
dev.off()

pdf("NSPCs_H3_expression.pdf")
get_hist_plot(my.npc.hist.H3)
dev.off()




