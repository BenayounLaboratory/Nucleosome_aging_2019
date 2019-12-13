setwd('/Volumes/BB_USC_1v2/Nucleosome_aging_project/Nucleosome_analyses_OLD_Mapping/Genome_Ontology/RANDOM/')
options(stringsAsFactors = F)

# 2019-11-04
# run parsed HOMER randomizations of detected nucleosomes (250 of median number of changed)
# 250 random samples of 5884 nucleosomes
my.cereb.gain.nucs  <- read.csv("2019-11-04_HOMER_GenOnt_Randomizations_Cerebellum_H3Nuc_GAIN.txt", sep = "\t", header = T)
my.cereb.lost.nucs  <- read.csv("2019-11-04_HOMER_GenOnt_Randomizations_Cerebellum_H3Nuc_LOST.txt", sep = "\t", header = T)
my.heart.gain.nucs  <- read.csv("2019-11-04_HOMER_GenOnt_Randomizations_Heart_H3Nuc_GAIN.txt", sep = "\t", header = T)
my.heart.lost.nucs  <- read.csv("2019-11-04_HOMER_GenOnt_Randomizations_Heart_H3Nuc_LOST.txt", sep = "\t", header = T)
my.liver.gain.nucs  <- read.csv("2019-11-04_HOMER_GenOnt_Randomizations_Liver_H3Nuc_GAIN.txt", sep = "\t", header = T)
my.liver.lost.nucs  <- read.csv("2019-11-04_HOMER_GenOnt_Randomizations_Liver_H3Nuc_LOST.txt", sep = "\t", header = T)
my.NPCs.gain.nucs   <- read.csv("2019-11-04_HOMER_GenOnt_Randomizations_NPCs_H3Nuc_GAIN.txt", sep = "\t", header = T)
my.NPCs.lost.nucs   <- read.csv("2019-11-04_HOMER_GenOnt_Randomizations_NPCs_H3Nuc_LOST.txt", sep = "\t", header = T)
my.OB.gain.nucs     <- read.csv("2019-11-04_HOMER_GenOnt_Randomizations_OlfactoryBulb_H3Nuc_GAIN.txt", sep = "\t", header = T)
my.OB.lost.nucs     <- read.csv("2019-11-04_HOMER_GenOnt_Randomizations_OlfactoryBulb_H3Nuc_LOST.txt", sep = "\t", header = T)

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

median(c(3348,5759,6557,7950,2458,4032,7777,6009,11178,1412))
# 5884



my.cereb.gain.nucs.enrich_Analysis  <- process_nucleosomes_stats("Cerebellum_GAIN", my.cereb.gain.nucs , 6557)
my.cereb.lost.nucs.enrich_Analysis  <- process_nucleosomes_stats("Cerebellum_LOSS", my.cereb.lost.nucs , 6009)
my.heart.gain.nucs.enrich_Analysis  <- process_nucleosomes_stats("Heart_GAIN", my.heart.gain.nucs , 3348)
my.heart.lost.nucs.enrich_Analysis  <- process_nucleosomes_stats("Heart_LOSS", my.heart.lost.nucs , 4032)
my.liver.gain.nucs.enrich_Analysis  <- process_nucleosomes_stats("Liver_GAIN", my.liver.gain.nucs , 5759)
my.liver.lost.nucs.enrich_Analysis  <- process_nucleosomes_stats("Liver_LOSS", my.liver.lost.nucs , 7777)
my.NPCs.gain.nucs.enrich_Analysis   <- process_nucleosomes_stats("NPCs_GAIN", my.NPCs.gain.nucs  , 2458)
my.NPCs.lost.nucs.enrich_Analysis   <- process_nucleosomes_stats("NPCs_LOSS", my.NPCs.lost.nucs  , 1412)
my.OB.gain.nucs.enrich_Analysis     <- process_nucleosomes_stats("OB_GAIN", my.OB.gain.nucs    , 7950)
my.OB.lost.nucs.enrich_Analysis     <- process_nucleosomes_stats("OB_LOSS", my.OB.lost.nucs    , 11178)

1/250 # 0.004, smallest pval possible

my.exclude <- c("centromeres", # NA ratio, exclude
                "gaps",        # NA ratio, exclude
                "rRNA",        # NA ratio, exclude
                "scRNA",       # NA ratio, exclude
                "snRNA",       # NA ratio, exclude
                "snoRNA",      # NA ratio, exclude
                "miRNA",       # 0 ratio, not meaningful
                "pseudo")      # 0 ratio, not meaningful

my.cereb.gain.nucs.enrich_Analysis.v2  <- my.cereb.gain.nucs.enrich_Analysis[  !(my.cereb.gain.nucs.enrich_Analysis$genomic_element %in% my.exclude),  ]
my.cereb.lost.nucs.enrich_Analysis.v2  <- my.cereb.lost.nucs.enrich_Analysis[  !(my.cereb.lost.nucs.enrich_Analysis$genomic_element %in% my.exclude),  ]
my.heart.gain.nucs.enrich_Analysis.v2  <- my.heart.gain.nucs.enrich_Analysis[  !(my.heart.gain.nucs.enrich_Analysis$genomic_element %in% my.exclude),  ]
my.heart.lost.nucs.enrich_Analysis.v2  <- my.heart.lost.nucs.enrich_Analysis[  !(my.heart.lost.nucs.enrich_Analysis$genomic_element %in% my.exclude),  ]
my.liver.gain.nucs.enrich_Analysis.v2  <- my.liver.gain.nucs.enrich_Analysis[  !(my.liver.gain.nucs.enrich_Analysis$genomic_element %in% my.exclude),  ]
my.liver.lost.nucs.enrich_Analysis.v2  <- my.liver.lost.nucs.enrich_Analysis[  !(my.liver.lost.nucs.enrich_Analysis$genomic_element %in% my.exclude),  ]
my.NPCs.gain.nucs.enrich_Analysis.v2   <- my.NPCs.gain.nucs.enrich_Analysis [  !(my.NPCs.gain.nucs.enrich_Analysis$genomic_element %in% my.exclude ),  ]
my.NPCs.lost.nucs.enrich_Analysis.v2   <- my.NPCs.lost.nucs.enrich_Analysis [  !(my.NPCs.lost.nucs.enrich_Analysis$genomic_element %in% my.exclude ),  ]
my.OB.gain.nucs.enrich_Analysis.v2     <- my.OB.gain.nucs.enrich_Analysis   [  !(my.OB.gain.nucs.enrich_Analysis$genomic_element %in% my.exclude   ),  ]
my.OB.lost.nucs.enrich_Analysis.v2     <- my.OB.lost.nucs.enrich_Analysis   [  !(my.OB.lost.nucs.enrich_Analysis$genomic_element %in% my.exclude   ),  ]

write.table(my.cereb.gain.nucs.enrich_Analysis.v2, file = paste0(Sys.Date(),"Cerebellum_Increased_H3_Nucleosomes", "logRatio_enrichment_Analysis.txt") , sep = "\t", quote = F)
write.table(my.cereb.lost.nucs.enrich_Analysis.v2, file = paste0(Sys.Date(),"Cerebellum_Decreased_H3_Nucleosomes", "logRatio_enrichment_Analysis.txt") , sep = "\t", quote = F)
write.table(my.heart.gain.nucs.enrich_Analysis.v2, file = paste0(Sys.Date(),"Heart_Increased_H3_Nucleosomes",      "logRatio_enrichment_Analysis.txt") , sep = "\t", quote = F)
write.table(my.heart.lost.nucs.enrich_Analysis.v2, file = paste0(Sys.Date(),"Heart_Decreased_H3_Nucleosomes",      "logRatio_enrichment_Analysis.txt") , sep = "\t", quote = F)
write.table(my.liver.gain.nucs.enrich_Analysis.v2, file = paste0(Sys.Date(),"Liver_Increased_H3_Nucleosomes",      "logRatio_enrichment_Analysis.txt") , sep = "\t", quote = F)
write.table(my.liver.lost.nucs.enrich_Analysis.v2, file = paste0(Sys.Date(),"Liver_Decreased_H3_Nucleosomes",      "logRatio_enrichment_Analysis.txt") , sep = "\t", quote = F)
write.table(my.NPCs.gain.nucs.enrich_Analysis.v2 , file = paste0(Sys.Date(),"NPCs_Increased_H3_Nucleosomes",       "logRatio_enrichment_Analysis.txt") , sep = "\t", quote = F)
write.table(my.NPCs.lost.nucs.enrich_Analysis.v2 , file = paste0(Sys.Date(),"NPCs_Decreased_H3_Nucleosomes",       "logRatio_enrichment_Analysis.txt") , sep = "\t", quote = F)
write.table(my.OB.gain.nucs.enrich_Analysis.v2   , file = paste0(Sys.Date(),"OB_Increased_H3_Nucleosomes",         "logRatio_enrichment_Analysis.txt") , sep = "\t", quote = F)
write.table(my.OB.lost.nucs.enrich_Analysis.v2   , file = paste0(Sys.Date(),"OB_Decreased_H3_Nucleosomes",         "logRatio_enrichment_Analysis.txt") , sep = "\t", quote = F)

plot_nucleosomes("Cerebellum Increased H3 Nucleosomes", my.cereb.gain.nucs.enrich_Analysis.v2 , my.gain = TRUE)
plot_nucleosomes("Cerebellum Decreased H3 Nucleosomes", my.cereb.lost.nucs.enrich_Analysis.v2 , my.gain = FALSE)
plot_nucleosomes("Heart Increased H3 Nucleosomes",      my.heart.gain.nucs.enrich_Analysis.v2 , my.gain = TRUE)
plot_nucleosomes("Heart Decreased H3 Nucleosomes",      my.heart.lost.nucs.enrich_Analysis.v2 , my.gain = FALSE)
plot_nucleosomes("Liver Increased H3 Nucleosomes",      my.liver.gain.nucs.enrich_Analysis.v2 , my.gain = TRUE)
plot_nucleosomes("Liver Decreased H3 Nucleosomes",      my.liver.lost.nucs.enrich_Analysis.v2 , my.gain = FALSE)
plot_nucleosomes("NPCs Increased H3 Nucleosomes",       my.NPCs.gain.nucs.enrich_Analysis.v2  , my.gain = TRUE)
plot_nucleosomes("NPCs Decreased H3 Nucleosomes",       my.NPCs.lost.nucs.enrich_Analysis.v2  , my.gain = FALSE)
plot_nucleosomes("OB Increased H3 Nucleosomes",         my.OB.gain.nucs.enrich_Analysis.v2    , my.gain = TRUE)
plot_nucleosomes("OB Decreased H3 Nucleosomes",         my.OB.lost.nucs.enrich_Analysis.v2    , my.gain = FALSE)


###############################################################
############              FUNCTIONS               #############
###############################################################

###############################################################
plot_nucleosomes <- function (my.title, my.enrich_Analysis.v2, my.gain = TRUE) {
  
  my.col = "#333399" ###loss
  
  if( my.gain ) {
    my.col = "#CC3333" 
  }
  
  pdf(paste0(Sys.Date(),my.title,"logRatio_enrichment.pdf"), height = 5, width = 5)
  par(oma=c(1,4,1,1))
  barplot(log2(my.enrich_Analysis.v2$Real_vs_simulated_ratio), 
          horiz = T, las = 2, col = my.col,
          names = my.enrich_Analysis.v2$genomic_element,
          main = my.title ,
          xlab = "Log2(Observed/Random Overlap)",
          xlim = c(-1.5,1.5) )
  abline(v=0)
  box()
  dev.off()
  
}

###############################################################
norm_to_nucs <- function(my.vector) {
  # divide by sampled nucleosomes
  return(my.vector/5884);
}


# test
# my.diff.nucs <- 6557
# my.nucs <- my.cereb.gain.nucs
# my.tissue <- "Cerebellum GAIN"

###############################################################
process_nucleosomes_stats <- function (my.tissue, my.nucs, my.diff.nucs) {
  
  # extract names, real vaues and simulated values
  my.genom.elem <- my.nucs[,1]
  my.real   <- my.nucs[,2]
  my.random <- my.nucs[,-c(1:2)]
  
  # normalize to number of nucleosomes to limit impact of set size
  # randoms all 5884 by design (median number of diff nucleosomes)
  # sum slightly over one since one nucleosome can overlap a boundary
  my.real.norm <- my.real/my.diff.nucs
  my.rand.norm <- apply(my.random, 2, norm_to_nucs)
  
  my.pval.enrich  <- rep(1,nrow(my.rand.norm))
  my.pval.deplete <- rep(1,nrow(my.rand.norm))
  my.ratio <- rep(1,nrow(my.rand.norm))
  
  pdf(paste0(Sys.Date(),"_ecdf_over_randoms_",my.tissue,".pdf"), height = 20, width = 15)
  par(oma=c(0.15,0.15,0.15,0.15))
  par(mfrow=c(5,4))
  
  # run over rows
  for (i in 1:nrow(my.rand.norm)) {
    
    my.ecdf <- ecdf(my.rand.norm[i,])

    my.pval.enrich[i]  <- 1-my.ecdf(my.real.norm[i])
    my.pval.deplete[i] <- my.ecdf(my.real.norm[i])
    
    my.ratio[i]  <- my.real.norm[i]/median(my.rand.norm[i,])
    
    
    hist(as.numeric(my.rand.norm[i,]), breaks = 50, xlim = c(0,1), 
         main = paste0(my.tissue, ", ", my.genom.elem[i]),
         xlab = "Fraction overlap with HOMER annotation")
    abline(v = my.real.norm[i], col = "red")

  }
  dev.off()
  
  my.enrich.analysis <- data.frame("genomic_element"    = my.genom.elem,
                                   "P_value_Enrichment" = my.pval.enrich,
                                   "P_value_Depletion"  = my.pval.deplete,
                                   "Real_vs_simulated_ratio" = my.ratio
                                   )
  return(my.enrich.analysis)
  
}
#############################################################################################################################################################################################
# log2(my.enrich.analysis$Real_vs_simulated_ratio)
# 
# plot(ecdf(my.rand.norm[6,]))
# 
# my.ecdf <- ecdf(my.cereb.gain.nucs[6,-c(1:2)])
# 
# my.ecdf(my.cereb.gain.nucs[6,2])
# 
# sd(as.numeric(my.cereb.gain.nucs[6,-c(1:2)]))

