setwd('/Volumes/BB_USC_1v2/Nucleosome_aging_project/Nucleosome_analyses_OLD_Mapping/GREAT/H3_nucleosomes')
options(stringsAsFactors=F)
library('pheatmap')
#source('parsing_functions.R')
source('parsing_functions_v2.R')


# 2016-08-31
# parse out GREAT results

# 2016-09-14
# get summaries

# 2016-11-17
# update using consensus changed nucleosomes

# 2017-04-07
# update output using increased contrast (log2/log10 enrichment?)

# 2019-06-07
# More filtering options

############
# my.cereb.gained <- read.csv("Cerebellum_DiNUP_DANPOS_GAINED.bed_GREAT.tsv", sep="\t",skip=4)
# my.cereb.gained <- my.cereb.gained[!is.na(my.cereb.gained$HyperBonfP),]
# sum(my.cereb.gained$HyperFdrQ < 0.05) # 1681
# #my.cereb.gained[my.cereb.gained$HyperFdrQ < 0.05,3]
# 
# my.cereb.loss <- read.csv("Cerebellum_DiNUP_DANPOS_LOST.bed_GREAT.tsv", sep="\t",skip=4)
# my.cereb.loss <- my.cereb.loss[!is.na(my.cereb.loss$HyperBonfP),]
# sum(my.cereb.loss$HyperFdrQ < 0.05) # 840
# #my.cereb.loss[my.cereb.loss$HyperFdrQ < 0.05,3]
# 
# 
# ############
# my.Heart.gained <- read.csv("Heart_DiNUP_DANPOS_GAINED.bed_GREAT.tsv", sep="\t",skip=4)
# my.Heart.gained <- my.Heart.gained[!is.na(my.Heart.gained$HyperBonfP),]
# sum(my.Heart.gained$HyperFdrQ < 0.05) # 1992
# #my.Heart.gained[my.Heart.gained$HyperFdrQ < 0.05,3]
# 
# my.Heart.loss <- read.csv("Heart_DiNUP_DANPOS_LOST.bed_GREAT.tsv", sep="\t",skip=4)
# my.Heart.loss <- my.Heart.loss[!is.na(my.Heart.loss$HyperBonfP),]
# sum(my.Heart.loss$HyperFdrQ < 0.05) # 2466
# #my.Heart.loss[my.Heart.loss$HyperFdrQ < 0.05,3]
# 
# 
# ############
# my.Liver.gained <- read.csv("Liver_DiNUP_DANPOS_GAINED.bed_GREAT.tsv", sep="\t",skip=4)
# my.Liver.gained <- my.Liver.gained[!is.na(my.Liver.gained$HyperBonfP),]
# sum(my.Liver.gained$HyperFdrQ < 0.05) # 2828
# #my.Liver.gained[my.Liver.gained$HyperFdrQ < 0.05,3]
# 
# my.Liver.loss <- read.csv("Liver_DiNUP_DANPOS_LOST.bed_GREAT.tsv", sep="\t",skip=4)
# my.Liver.loss <- my.Liver.loss[!is.na(my.Liver.loss$HyperBonfP),]
# sum(my.Liver.loss$HyperFdrQ < 0.05) # 2847
# #my.Liver.loss[my.Liver.loss$HyperFdrQ < 0.05,3]
# 
# 
# ############
# my.OB.gained <- read.csv("OlfactoryBulb_DiNUP_DANPOS_GAINED.bed_GREAT.tsv", sep="\t",skip=4)
# my.OB.gained <- my.OB.gained[!is.na(my.OB.gained$HyperBonfP),]
# sum(my.OB.gained$HyperFdrQ < 0.05) # 1752
# #my.OB.gained[my.OB.gained$HyperFdrQ < 0.05,3]
# 
# my.OB.loss <- read.csv("OlfactoryBulb_DiNUP_DANPOS_LOST.bed_GREAT.tsv", sep="\t",skip=4)
# my.OB.loss <- my.OB.loss[!is.na(my.OB.loss$HyperBonfP),]
# sum(my.OB.loss$HyperFdrQ < 0.05) # 1655
# #my.OB.loss[my.OB.loss$HyperFdrQ < 0.05,3]
# 
# 
# ############
# my.NPCs.gained <- read.csv("NPCs_DiNUP_DANPOS_GAINED.bed_GREAT.tsv", sep="\t",skip=4)
# my.NPCs.gained <- my.NPCs.gained[!is.na(my.NPCs.gained$HyperBonfP),]
# sum(my.NPCs.gained$HyperFdrQ < 0.05) # 1760
# #my.NPCs.gained[my.NPCs.gained$HyperFdrQ < 0.05,3]
# 
# my.NPCs.loss <- read.csv("NPCs_DiNUP_DANPOS_LOST.bed_GREAT.tsv", sep="\t",skip=4)
# my.NPCs.loss <- my.NPCs.loss[!is.na(my.NPCs.loss$HyperBonfP),]
# sum(my.NPCs.loss$HyperFdrQ < 0.05) # 1226
# #my.NPCs.loss[my.NPCs.loss$HyperFdrQ < 0.05,3]
# 
# 
# 
# ####################################################################################################################################################
# # filter relevant columns
# my.cereb.gained.f <- my.cereb.gained[,c("X..Ontology","ID","Desc","RegionFoldEnrich", "HyperP","HyperFdrQ")]
# my.cereb.loss.f <- my.cereb.loss[,c("X..Ontology","ID","Desc","RegionFoldEnrich", "HyperP","HyperFdrQ")]
# my.Heart.gained.f <- my.Heart.gained[,c("X..Ontology","ID","Desc","RegionFoldEnrich", "HyperP","HyperFdrQ")]
# my.Heart.loss.f <- my.Heart.loss[,c("X..Ontology","ID","Desc","RegionFoldEnrich", "HyperP","HyperFdrQ")]
# my.Liver.gained.f <- my.Liver.gained[,c("X..Ontology","ID","Desc","RegionFoldEnrich", "HyperP","HyperFdrQ")]
# my.Liver.loss.f <- my.Liver.loss[,c("X..Ontology","ID","Desc","RegionFoldEnrich", "HyperP","HyperFdrQ")]
# my.OB.gained.f <- my.OB.gained[,c("X..Ontology","ID","Desc","RegionFoldEnrich", "HyperP","HyperFdrQ")]
# my.OB.loss.f <- my.OB.loss[,c("X..Ontology","ID","Desc","RegionFoldEnrich", "HyperP","HyperFdrQ")]
# my.NPCs.gained.f <- my.NPCs.gained[,c("X..Ontology","ID","Desc","RegionFoldEnrich", "HyperP","HyperFdrQ")]
# my.NPCs.loss.f <- my.NPCs.loss[,c("X..Ontology","ID","Desc","RegionFoldEnrich", "HyperP","HyperFdrQ")]
# 
# 
# # rename columns for merging
# colnames(my.Heart.gained.f) <- c("Ontology","ID","Description","Heart_RegionFoldEnrich", "Heart_HyperP","Heart_HyperFdrQ")
# colnames(my.Liver.gained.f) <- c("Ontology","ID","Description","Liver_RegionFoldEnrich", "Liver_HyperP","Liver_HyperFdrQ")
# colnames(my.cereb.gained.f) <- c("Ontology","ID","Description","Cerebellum_RegionFoldEnrich", "Cerebellum_HyperP","Cerebellum_HyperFdrQ")
# colnames(my.OB.gained.f) <- c("Ontology","ID","Description","OB_RegionFoldEnrich", "OB_HyperP","OB_HyperFdrQ")
# colnames(my.NPCs.gained.f) <- c("Ontology","ID","Description","NSPCs_RegionFoldEnrich", "NSPCs_HyperP","NSPCs_HyperFdrQ")
# 
# colnames(my.Heart.loss.f) <- c("Ontology","ID","Description","Heart_RegionFoldEnrich", "Heart_HyperP","Heart_HyperFdrQ")
# colnames(my.Liver.loss.f) <- c("Ontology","ID","Description","Liver_RegionFoldEnrich", "Liver_HyperP","Liver_HyperFdrQ")
# colnames(my.cereb.loss.f) <-c("Ontology","ID","Description","Cerebellum_RegionFoldEnrich", "Cerebellum_HyperP","Cerebellum_HyperFdrQ")
# colnames(my.OB.loss.f) <- c("Ontology","ID","Description","OB_RegionFoldEnrich", "OB_HyperP","OB_HyperFdrQ")
# colnames(my.NPCs.loss.f) <- c("Ontology","ID","Description","NSPCs_RegionFoldEnrich", "NSPCs_HyperP","NSPCs_HyperFdrQ")
# 
# 
# ##########################################################################
# ################   merge gained nucleosomes   ############################
# ##########################################################################
# 
# my.merge.gain <- merge(my.Heart.gained.f,my.Liver.gained.f, by='ID', all = T)
# my.merge.gain2 <- merge(my.cereb.gained.f,my.OB.gained.f, by='ID', all = T)
# my.merge.gain3 <- merge(my.merge.gain,my.merge.gain2, by='ID', all = T)
# my.merge.gain4 <- merge(my.merge.gain3,my.NPCs.gained.f, by='ID', all = T)
# dim(my.merge.gain4)
# # 23153    26
# 
# 
# # compile clean result matrix
# my.results <- data.frame(matrix(0,dim(my.merge.gain4)[1],18))
# colnames(my.results) <- c("Ontology","ID","Description","Heart_RegionFoldEnrich", "Heart_HyperP","Heart_HyperFdrQ",
#                           "Liver_RegionFoldEnrich", "Liver_HyperP","Liver_HyperFdrQ",
#                           "Cerebellum_RegionFoldEnrich", "Cerebellum_HyperP","Cerebellum_HyperFdrQ",
#                           "OB_RegionFoldEnrich", "OB_HyperP","OB_HyperFdrQ",
#                           "NSPCs_RegionFoldEnrich", "NSPCs_HyperP","NSPCs_HyperFdrQ")
# 
# my.descs <- grep("Description",colnames(my.merge.gain4))
# my.onto <- grep("Ontology",colnames(my.merge.gain4))
# 
# for (i in 1:dim(my.merge.gain4)[1]) {
#   
#   # get descriptions
#   my.ok <- which(!is.na( my.merge.gain4[i,my.descs]))
#   my.results$Description[i] <- my.merge.gain4[i,my.descs[my.ok[1]]]
#   my.results$Ontology[i] <- my.merge.gain4[i,my.onto[my.ok[1]]]
#   my.results$ID[i] <- my.merge.gain4[i,"ID"]
#   
#   my.results[i,4:18] <- my.merge.gain4[i,-c(1, my.descs,my.onto)]
#   
#   # replace NA by 0 enrich, and 1 p-value
#   my.results[i,c(4,7,10,13,16)[ which(is.na(my.results[i,c(4,7,10,13,16)] ) ) ] ]<- 0
#   my.results[i,c(5:6,8:9,11:12,14:15,17:18) [which(is.na(my.results[i,c(5:6,8:9,11:12,14:15,17:18)]))] ] <- 1
#   
#   
# }
# 
# # extract those with at least one tissue at FDR 0.05
# my.sigs.nums <- apply(my.results[,c("Heart_HyperFdrQ","Liver_HyperFdrQ","Cerebellum_HyperFdrQ","OB_HyperFdrQ","NSPCs_HyperFdrQ")] < 0.05,1, sum)
# sum(my.sigs.nums >= 1) # 5295
# 
# write.table(my.results[my.sigs.nums >= 1,], sep="\t",row.names=F,quote=F, 
#             file = paste(Sys.Date(),"H3_gained_nucleosomes_GREAT_significant_in_1ormore_tissues_FDR_0.05.txt"))
# 
# 
# my.sigs.nums2 <- apply(my.results[,c("Heart_HyperFdrQ","Liver_HyperFdrQ","Cerebellum_HyperFdrQ","OB_HyperFdrQ","NSPCs_HyperFdrQ")] < 0.01,1, sum)
# sum(my.sigs.nums2 >= 1) # 3750
# 
# write.table(my.results[my.sigs.nums >= 1,], sep="\t",row.names=F,quote=F, 
#             file = paste(Sys.Date(),"H3_gained_nucleosomes_GREAT_significant_in_1ormore_tissues_FDR_0.01.txt"))
# 
# 
# # my.enrich.res <- c("Heart_HyperFdrQ","Liver_HyperFdrQ","Cerebellum_HyperFdrQ","OB_HyperFdrQ","NSPCs_HyperFdrQ")
# 
# ### restrict to FDR 5 and less
# my.res2 <- my.results[my.sigs.nums >= 1,]
# save(my.res2, file = "2016-11-17_summary_matrix_consensus_GAINS_H3_nucleosomes.RData")

###########################
# load from file
load('2016-11-17_summary_matrix_consensus_GAINS_H3_nucleosomes.RData')
my.enrich.res <- c("Heart_RegionFoldEnrich","Liver_RegionFoldEnrich","Cerebellum_RegionFoldEnrich","OB_RegionFoldEnrich","NSPCs_RegionFoldEnrich")
  
#### 1. MGI Expression: Detected
plot_gains_or_losses(my.res2, "MGI Expression: Detected", 3, c(0.05,1e-3,1e-6), "GAINED", "H3_nucleosomes")

#### 2. MSigDB miRNA Motifs
plot_gains_or_losses(my.res2, "MSigDB miRNA Motifs",3, c(0.05,1e-3,1e-6), "GAINED", "H3_nucleosomes")

#### 3. MSigDB Predicted Promoter Motifs
plot_gains_or_losses(my.res2, "MSigDB Predicted Promoter Motifs", 3, c(0.05,1e-3,1e-6), "GAINED", "H3_nucleosomes")

#### 4. MSigDB Perturbation
plot_gains_or_losses(my.res2, "MSigDB Perturbation", 3, c(0.05,1e-3,1e-6), "GAINED", "H3_nucleosomes")

#### 5. Disease Ontology
plot_gains_or_losses(my.res2, "Disease Ontology", 3, c(0.05,1e-3,1e-6), "GAINED", "H3_nucleosomes")

#### 6. GO Biological Process
plot_gains_or_losses(my.res2, "GO Biological Process", 3, c(0.05,1e-3,1e-6), "GAINED", "H3_nucleosomes")

#### 7. GO Molecular Function
plot_gains_or_losses(my.res2, "GO Molecular Function", 3, c(0.05,1e-3,1e-6), "GAINED", "H3_nucleosomes")

#### 8. GO Cellular Component
plot_gains_or_losses(my.res2, "GO Cellular Component", 3, c(0.05,1e-3,1e-6), "GAINED", "H3_nucleosomes")

#### 9. Mouse Phenotype
plot_gains_or_losses(my.res2, "Mouse Phenotype", 3, c(0.05,1e-3,1e-6), "GAINED", "H3_nucleosomes")


##########################################################################
################    merge Lost nucleosomes    ############################
##########################################################################

# my.merge.loss <- merge(my.Heart.loss.f,my.Liver.loss.f, by='ID', all = T)
# my.merge.loss2 <- merge(my.cereb.loss.f,my.OB.loss.f, by='ID', all = T)
# my.merge.loss3 <- merge(my.merge.loss,my.merge.loss2, by='ID', all = T)
# my.merge.loss4 <- merge(my.merge.loss3,my.NPCs.loss.f, by='ID', all = T)
# dim(my.merge.loss4)
# #[1] 24113    26
# 
# # compile clean result matrix
# my.results <- data.frame(matrix(0,dim(my.merge.loss4)[1],18))
# colnames(my.results) <- c("Ontology","ID","Description","Heart_RegionFoldEnrich", "Heart_HyperP","Heart_HyperFdrQ",
#                           "Liver_RegionFoldEnrich", "Liver_HyperP","Liver_HyperFdrQ",
#                           "Cerebellum_RegionFoldEnrich", "Cerebellum_HyperP","Cerebellum_HyperFdrQ",
#                           "OB_RegionFoldEnrich", "OB_HyperP","OB_HyperFdrQ",
#                           "NSPCs_RegionFoldEnrich", "NSPCs_HyperP","NSPCs_HyperFdrQ")
# 
# my.descs <- grep("Description",colnames(my.merge.loss4))
# my.onto <- grep("Ontology",colnames(my.merge.loss4))
# 
# for (i in 1:dim(my.merge.loss4)[1]) {
#   
#   # get descriptions
#   my.ok <- which(!is.na( my.merge.loss4[i,my.descs]))
#   my.results$Description[i] <- my.merge.loss4[i,my.descs[my.ok[1]]]
#   my.results$Ontology[i] <- my.merge.loss4[i,my.onto[my.ok[1]]]
#   my.results$ID[i] <- my.merge.loss4[i,"ID"]
#   
#   my.results[i,4:18] <- my.merge.loss4[i,-c(1, my.descs,my.onto)]
#   
#   # replace NA by 0 enrich, and 1 p-value
#   my.results[i,c(4,7,10,13,16)[ which(is.na(my.results[i,c(4,7,10,13,16)] ) ) ] ]<- 0
#   my.results[i,c(5:6,8:9,11:12,14:15,17:18) [which(is.na(my.results[i,c(5:6,8:9,11:12,14:15,17:18)]))] ] <- 1
#   
#   
# }
# 
# # extract those with at least one tissue at FDR 0.05
# my.sigs.nums <- apply(my.results[,c("Heart_HyperFdrQ","Liver_HyperFdrQ","Cerebellum_HyperFdrQ","OB_HyperFdrQ","NSPCs_HyperFdrQ")] < 0.05,1, sum)
# sum(my.sigs.nums >= 1) # 5184
# 
# write.table(my.results[my.sigs.nums >= 1,], sep="\t",row.names=F,quote=F, 
#             file = paste(Sys.Date(),"H3_lost_nucleosomes_GREAT_significant_in_1ormore_tissues_FDR_0.05.txt"))
# 
# 
# my.sigs.nums2 <- apply(my.results[,c("Heart_HyperFdrQ","Liver_HyperFdrQ","Cerebellum_HyperFdrQ","OB_HyperFdrQ","NSPCs_HyperFdrQ")] < 0.01,1, sum)
# sum(my.sigs.nums2 >= 1) # 3706
# 
# write.table(my.results[my.sigs.nums2 >= 1,], sep="\t",row.names=F,quote=F, 
#             file = paste(Sys.Date(),"H3_lost_nucleosomes_GREAT_significant_in_1ormore_tissues_FDR_0.01.txt"))
# 
# ### restrict to FDR 5 and less
# my.res2 <- my.results[my.sigs.nums >= 1,]
# save(my.res2, file = "2016-11-17_summary_matrix_consensus_LOSSES_H3_nucleosomes.RData")


###########################
load('2016-11-17_summary_matrix_consensus_LOSSES_H3_nucleosomes.RData')

#### 1. MGI Expression: Detected
plot_gains_or_losses(my.res2, "MGI Expression: Detected", 3, c(0.05,1e-3,1e-6), "LOST", "H3_nucleosomes")

#### 2. MSigDB miRNA Motifs
plot_gains_or_losses(my.res2, "MSigDB miRNA Motifs", 3, c(0.05,1e-3,1e-6), "LOST", "H3_nucleosomes")

#### 3. MSigDB Predicted Promoter Motifs
plot_gains_or_losses(my.res2, "MSigDB Predicted Promoter Motifs", 3, c(0.05,1e-3,1e-6), "LOST", "H3_nucleosomes")

#### 4. MSigDB Perturbation
plot_gains_or_losses(my.res2, "MSigDB Perturbation", 3, c(0.05,1e-3,1e-6), "LOST", "H3_nucleosomes")

#### 5. Disease Ontology
plot_gains_or_losses(my.res2, "Disease Ontology", 3, c(0.05,1e-3,1e-6), "LOST", "H3_nucleosomes")

#### 6. GO Biological Process
plot_gains_or_losses(my.res2, "GO Biological Process", 3, c(0.05,1e-3,1e-6), "LOST", "H3_nucleosomes")

#### 7. GO Molecular Function
plot_gains_or_losses(my.res2, "GO Molecular Function", 3, c(0.05,1e-3,1e-6), "LOST", "H3_nucleosomes")

#### 8. GO Cellular Component
plot_gains_or_losses(my.res2, "GO Cellular Component", 3, c(0.05,1e-3,1e-6), "LOST", "H3_nucleosomes")

#### 9. Mouse Phenotype
plot_gains_or_losses(my.res2, "Mouse Phenotype", 3, c(0.05,1e-3,1e-6), "LOST", "H3_nucleosomes")
