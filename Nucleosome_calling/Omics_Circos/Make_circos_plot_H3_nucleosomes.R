setwd('/Volumes/BB_USC_1/Nucleosome_aging_project/Nucleosome_analyses_OLD_Mapping/Circos/')
options(stringsAsFactors = FALSE);
library("OmicCircos")

source('omics_circos_functions_v2.R')

# 2019-02-07
# plot changing nucleosomes

###########
## construct mm9 scaffold, letting an opening so the legend can be written in
my.mm9 <- read.table('~/Softwares/bedtools2.26/genomes/mouse.mm9.genome', sep="\t")
my.mm9 <- data.frame('chr' = my.mm9$V1,
                     'start'= rep(0,length(my.mm9$V1)),
                     'end'= my.mm9$V2,
                     'bla1'= rep(NA,length(my.mm9$V1)),
                     'bla2'= rep(NA,length(my.mm9$V1))
)
my.mm9 <- my.mm9[1:21,]

mm9.db <- segAnglePo(my.mm9, seg=my.mm9$chr, angle.end = 350);


### process datasets
# Liver
my.h3.liver <- read.table('2017-03-20_DiNuP_DANPOS_Liver_H3_aging.bed',sep="\t",header=F)
my.h3.liver.in <- make_omics_Nucs(my.h3.liver )        

# Heart
my.h3.heart <- read.table('2017-03-20_DiNuP_DANPOS_Heart_H3_aging.bed',sep="\t",header=F)
my.h3.heart.in <- make_omics_Nucs(my.h3.heart )

# Cerebellum
my.h3.cereb <- read.table('2017-03-20_DiNuP_DANPOS_Cerebellum_H3_aging.bed',sep="\t",header=F)
my.h3.cereb.in <- make_omics_Nucs(my.h3.cereb )

# OB
my.h3.ob    <- read.table('2017-03-20_DiNuP_DANPOS_OB_H3_aging.bed',sep="\t",header=F)
my.h3.ob.in    <- make_omics_Nucs(my.h3.ob    )     

# NPCs
my.h3.npcs  <- read.table('2017-03-20_DiNuP_DANPOS_NPCs_H3_aging.bed',sep="\t",header=F)
my.h3.npcs.in  <- make_omics_Nucs(my.h3.npcs  )


######################################################################################################

##### generate plots

my.nuc.colors <- c("firebrick1","deepskyblue2")

pdf(paste(Sys.Date(),"Circos_all_DE_H3_aging.pdf", sep="_"), width = 8, height = 8)
par(mar=c(0, 0, 0, 0));

plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="");
circos(R=410, cir=mm9.db, W=15, type="chr", print.chr.lab=TRUE, scale=TRUE, col = "black");

# Heart
circos(R=360, cir=mm9.db, W=50,  mapping=my.h3.heart.in,
       col.v=3, type="b2", B=T, lwd=0.5, cutoff=0, col=my.nuc.colors);

# Liver
circos(R=300, cir=mm9.db, W=50,  mapping=my.h3.liver.in,
       col.v=3, type="b2", B=T, lwd=0.5, cutoff=0, col=my.nuc.colors);

# Cerebellum
circos(R=240, cir=mm9.db, W=50,  mapping=my.h3.cereb.in,
       col.v=3, type="b2", B=T, lwd=0.5, cutoff=0, col=my.nuc.colors);

# OB
circos(R=190, cir=mm9.db, W=50,  mapping=my.h3.ob.in,
       col.v=3, type="b2", B=T, lwd=0.5, cutoff=0, col=my.nuc.colors);

# NPCs
circos(R=140, cir=mm9.db, W=50,  mapping=my.h3.npcs.in,
       col.v=3, type="b2", B=T, lwd=0.5, cutoff=0, col=my.nuc.colors);
dev.off()
