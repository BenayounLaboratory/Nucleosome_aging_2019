# 2017-04-07
# update output using increased contrast (log2/log10 enrichment?)
# was plotting pvalues => plot log2 enrich fold
# plot_gains(my.res2, "MGI Expression: Detected", 1, c(0.05,1e-3,1e-6))
# my.num.sigs <- 1
# fdr.thrs <- c(0.05,1e-3,1e-6)
# my.ontology <- "MGI Expression: Detected"

my.enrich.res <- c("Heart_RegionFoldEnrich","Liver_RegionFoldEnrich","Cerebellum_RegionFoldEnrich","OB_RegionFoldEnrich","NSPCs_RegionFoldEnrich")


plot_gains_or_losses <- function(my.res2, my.ontology, my.num.sigs, fdr.thrs, my.direction, my.mark, my.enrich.results = my.enrich.res) {
  
  for ( i in 1:length(fdr.thrs)) {
    
    my.sigs.nums <- apply(my.res2[,c("Heart_HyperFdrQ","Liver_HyperFdrQ","Cerebellum_HyperFdrQ","OB_HyperFdrQ","NSPCs_HyperFdrQ")] < fdr.thrs[i] ,1, sum)
    
    my.sig.gain1 <- my.res2[my.sigs.nums >= my.num.sigs,]
    
    #my.ontology:  "MGI Expression: Detected"
    my.mgi.exp <- my.sig.gain1[my.sig.gain1$Ontology %in% my.ontology,]
    rownames(my.mgi.exp) <- paste(my.mgi.exp$ID,my.mgi.exp$Description,1:length(my.mgi.exp$ID),sep="_")
    
    print(paste("FDR",fdr.thrs[i],":",dim(my.mgi.exp)[1]))
    
    my.res.mat <- log2(my.mgi.exp[,my.enrich.results] + 0.1) # log2 fold enrichments
    
    my.max <- max(abs(my.res.mat))
    
    #     my.scaling <- rbind(rep(my.max, 5),
    #                         rep(-my.max, 5))
    #     rownames(my.scaling) <- c("MAX","MIN")
    #     colnames(my.scaling) <- my.enrich.results
    
    paletteLength <- 200
    
    my.heat.colors <- colorRampPalette(c("gray88","gray92","gray95","white","firebrick1","firebrick2","firebrick4"))(paletteLength)
    if ( (tolower(my.direction) %in% 'lost') || (tolower(my.direction) %in% 'decreased') ) {
      my.heat.colors <- colorRampPalette(c("gray88","gray92","gray95","white","deepskyblue","dodgerblue2","dodgerblue4"))(paletteLength)
    }
    
    # see answer:http://stackoverflow.com/questions/31677923/set-0-point-for-pheatmap-in-r
    # length(breaks) == length(paletteLength) + 1
    # use floor and ceiling to deal with even/odd length pallettelengths
    myBreaks <- c(seq(min(my.res.mat), 0, length.out=ceiling(paletteLength/2) + 1), seq(max(my.res.mat)/paletteLength, max(my.res.mat), length.out=floor(paletteLength/2)))
    
    
    my.pdfname <-paste(Sys.Date(),my.direction,my.mark,"RegionFoldEnrich_heatmap",my.ontology,"significant_in", my.num.sigs, "or_more_FDR",fdr.thrs[i],".pdf", sep="_")
    pdf(my.pdfname, onefile=F, height = max(length(my.mgi.exp$ID)/9, 6), width=15)
    pheatmap(my.res.mat,cluster_cols = F,cluster_row = T,
             color = my.heat.colors, breaks = myBreaks,
             cellwidth = 40, cellheight = 5,
             fontsize_row = 5, border_color = NA)
    dev.off()
    
    my.txt.name <-paste(Sys.Date(),my.direction,my.mark,"RegionFoldEnrich_Table",my.ontology,"significant_in", my.num.sigs, "or_more_FDR",fdr.thrs[i],".txt", sep="_")
    
    write.table(my.res.mat,file = my.txt.name, quote = F, sep = "\t")
    
  }
  
  
}
