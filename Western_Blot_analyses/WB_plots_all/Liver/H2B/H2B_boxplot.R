
#############################################################################################
setwd('/Users/Yilin/Desktop/H3_resubmit_WB/WB_plots_all/Liver/H2B')

my.data <- read.table('H2B_liver_norm5m_median_5+21m.txt', sep = "\t", header = T)

my.colnames <- c("5m", "21m")
colnames(my.data)<- my.colnames
pdf(paste(Sys.Date(),"H2B_liver_norm5m_median_5+21m.pdf", sep = "_"),height = 6, width = 3.5)
boxplot(my.data, ylim = c(0,3), 
        col = c( "darkslategray1", "dodgerblue2"),
        ylab = "Relative Protein Intesity Normalized to Vinculin",main="H2B with Aging Liver", outline = F)
points(jitter(rep(1,10), factor = 4),my.data$`5m`, pch = 16, col = c(rep("magenta4",5),rep("mediumorchid2",5)), cex = 2)
points(jitter(rep(2,10), factor = 4),my.data$`21m`, pch = 16, col = c(rep("magenta4",5),rep("mediumorchid2",5)), cex = 2)
#points(jitter(rep(2,13), factor = 3),my.data$`21m`, pch = 16, col = c(rep("red",5),rep("blue",5),rep("black",3)))
#points(jitter(rep(3,13), factor = 3),my.data$`12m`, pch = 16, col = c(rep("black",2)))
#points(jitter(rep(4,13), factor = 3),my.data$`25m`, pch = 16, col = c(rep("black",2)))
legend("topright",inset = 0.02,legend=c("cohort 1", "cohort 2"),fill = c("magenta4","mediumorchid2"))
dev.off()

wilcox.test(my.data$`5m`,my.data$`21m`)