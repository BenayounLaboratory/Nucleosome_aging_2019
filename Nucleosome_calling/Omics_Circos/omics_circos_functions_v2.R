options(stringsAsFactors=F)


# extract coordinates from DiffBind/DESeq outputs
get_coords <- function(DEseq.res){
  coords.raw <- strsplit(rownames(DEseq.res), '-', fixed = TRUE)
  my.coords <- data.frame(matrix("",length(coords.raw),3))
  for (i in 1:length(coords.raw)){
    my.coords[i,] <- coords.raw[[i]]
  }
  colnames(my.coords) <- c('chr','start','end')
  
  return(my.coords)
}


# generate an input dataframe for omics circos for only significant regions
make_omics_fromDE_chrom <- function(my.height.process, p_thres = 0.1) {
  
  my.height.process.sigs <- data.frame(my.height.process[my.height.process$padj < p_thres,])
  
  my.coords <- get_coords(my.height.process.sigs)
  
  po <- 0.5*(as.numeric(my.coords$start)+as.numeric(my.coords$end))
  my.chrom.aging.in <- data.frame(cbind(my.coords$chr,
                                        po,
                                        data.frame(my.height.process.sigs[,c(7,6,2)])))
  colnames(my.chrom.aging.in)[1] <- 'chr'
  
  return(my.chrom.aging.in)
}

# generate an input dataframe for omics circos for only significant regions
make_omics_Nucs <- function(my.h3) {
    my.chrom.aging.in <- data.frame('chr' = my.h3$V1,
                                  'po' = 0.5*(my.h3$V2 + my.h3$V3),
                                  'logFC' = log2(my.h3$V5))
  
  
  return(my.chrom.aging.in)
}


###########
# my.breadth.process <- my.liver.breadth.process
# my.breadth.homer <- '/Volumes/MyBook_3/BD_aging_project/SMITE_data_intergration/ChIP-seq_data_files/Merged_ALL_AGES_MERGED_Liver_H3K4me3.PARSED_INTERSECTIONS.xls'
# p_thres = 0.05

make_omics_fromDE_chrom_breadth <- function(my.breadth.process, my.breadth.homer, p_thres = 0.1) {
  
  my.breadth.process.sigs <- data.frame(my.breadth.process[my.breadth.process$padj < p_thres,])
  my.breadth.process.sigs$PeakName <- rownames(my.breadth.process.sigs)
  my.breadths.homer.data <- read.csv(my.breadth.homer, header = T, sep = "\t")
  
  my.merged <- merge(my.breadth.process.sigs,my.breadths.homer.data[,c(1:4)], by.x = 'PeakName', by.y = 'Peak_Name')

  my.coords <- my.merged[,c('Chr','Start','End')]
  colnames(my.coords) <- c('chr','start','end')
  
  po <- 0.5*(as.numeric(my.coords$start)+as.numeric(my.coords$end))
  my.chrom.aging.in <- data.frame(cbind(my.coords$chr,
                                        po,
                                        data.frame(my.breadth.process.sigs[,c(7,6,2)])))
  colnames(my.chrom.aging.in)[1] <- 'chr'
  
  return(my.chrom.aging.in)
}

# generate an input dataframe for omics circos for only significant regions
make_omics_Nucs <- function(my.h3) {
  my.chrom.aging.in <- data.frame('chr' = my.h3$V1,
                                  'po' = 0.5*(my.h3$V2 + my.h3$V3),
                                  'logFC' = log2(my.h3$V5))
  
  
  return(my.chrom.aging.in)
}
# generate an input dataframe for omics circos for only significant genes (FDR 0.1)
# my.RNA.process <- my.liver.RNAseq.process
# p_thres = 0.05
# genes.coords = my.gene.coords

make_omics_fromDE_RNA <- function(my.RNA.process, genes.coords = my.gene.coords, p_thres = 0.05) {
  
  my.RNA.process.sigs <- data.frame(my.RNA.process[[2]][my.RNA.process[[1]]$padj < p_thres,])
  my.RNA.process.sigs$GeneName <- rownames(my.RNA.process.sigs)
  
  my.RNA.with.coords <- merge(my.RNA.process.sigs, my.gene.coords[,c("Chr","Start","End","Gene.Name")], by.x = 'GeneName', by.y = "Gene.Name" )
  
  po <- 0.5*(as.numeric(my.RNA.with.coords$Start)+as.numeric(my.RNA.with.coords$End))
  
  my.non.data.cols <- colnames(my.RNA.with.coords) %in% c("Chr","Start","End","GeneName")
  
  my.RNA.aging.in <- data.frame(cbind(my.RNA.with.coords$Chr,
                                        po,
                                        data.frame(my.RNA.with.coords[,!my.non.data.cols])))
  colnames(my.RNA.aging.in)[1] <- 'chr'
  
  return(my.RNA.aging.in)
}

# based on variable genes
make_omics_fromDE_RNA_v2 <- function(my.RNA.process, genes.coords = my.gene.coords, gene_vars = 1000) {
  
  # get absolute values of logFC
  my.logFCs <- abs(my.RNA.process[[1]]$log2FoldChange)
  my.sort.order <- sort(my.logFCs, decreasing = T, index.return = T)
  
  my.RNA.process.sigs <- data.frame(my.RNA.process[[2]][my.sort.order$ix[1:gene_vars],])
  my.RNA.process.sigs$GeneName <- rownames(my.RNA.process.sigs)
  
  my.RNA.with.coords <- merge(my.RNA.process.sigs, my.gene.coords[,c("Chr","Start","End","Gene.Name")], by.x = 'GeneName', by.y = "Gene.Name" )
  
  po <- 0.5*(as.numeric(my.RNA.with.coords$Start)+as.numeric(my.RNA.with.coords$End))
  
  my.non.data.cols <- colnames(my.RNA.with.coords) %in% c("Chr","Start","End","GeneName")
  
  my.RNA.aging.in <- data.frame(cbind(my.RNA.with.coords$Chr,
                                      po,
                                      data.frame(my.RNA.with.coords[,!my.non.data.cols])))
  colnames(my.RNA.aging.in)[1] <- 'chr'
  
  return(my.RNA.aging.in)
}
