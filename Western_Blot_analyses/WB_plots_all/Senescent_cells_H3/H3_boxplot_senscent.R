
#############################################################################################
setwd('/Users/Yilin/Desktop/WB_plots/Senescent_cells_H3')

my.data <- read.table('H3_not_norm.txt', sep = "\t", header = T)

my.colnames <- c("Proliferating", "Senescent")
colnames(my.data)<- my.colnames
pdf(paste(Sys.Date(),"H3_not_norm.pdf", sep = "_"),height = 6, width = 3.5)
boxplot(my.data, ylim = c(0,2), 
        col = c( "darkslategray1", "dodgerblue2"),
        ylab = "Relative Protein Intesity Normalized to Vinculin",main="H3_17kD with Cells Normalize to Control Average", outline = F)
points(jitter(rep(1,7), factor = 3),my.data$`Proliferating`, pch = 16, col = c(rep("magenta4",3),rep("mediumorchid2",3),,rep("darkviolet",1)), cex = 2)
points(jitter(rep(2,7), factor = 3),my.data$`Senescent`, pch = 16, col = c(rep("magenta4",3),rep("mediumorchid2",3),,rep("darkviolet",1)), cex = 2)
#points(jitter(rep(3,13), factor = 3),my.data$`12m`, pch = 16, col = c(rep("maroon1",1)),cex = 2)
#points(jitter(rep(4,13), factor = 3),my.data$`25m`, pch = 16, col = c(rep("black",2)))
legend("topright",inset = 0.02,legend=c("IMR90", "WI38","Human Fibroblasts"),fill = c("magenta4","mediumorchid2","darkviolet"))
dev.off()

wilcox.test(my.data$`Proliferating`,my.data$`Senescent`)

#wilcox.test(my.data$`5m`,my.data$`12m`)

#wilcox.test(my.data$`12m`,my.data$`21m`)

#wilcox.test(my.data$`12m`,my.data$`25m`)

#wilcox.test(my.data$`5m`,my.data$`25m`)

#wilcox.test(my.data$`21m`,my.data$`25m`)

