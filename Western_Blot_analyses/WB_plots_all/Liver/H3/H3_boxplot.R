
#############################################################################################
setwd('/Users/Yilin/Desktop/Uncropped_WB_pics/Liver/H3')


my.data <- read.table('H3_liver_norm5m_median_5+21m.txt', sep = "\t", header = T)

my.colnames <- c("5m", "21m")
colnames(my.data)<- my.colnames
pdf(paste(Sys.Date(),"H3_liver_norm5m_median_5+21m.pdf", sep = "_"),height = 6, width = 6)
boxplot(my.data, ylim = c(0,4), 
        col = c( "darkslategray1", "dodgerblue2"),
        ylab = "Relative Protein Intesity Normalized to Vinculin",main="H3_17kD with Aging Liver Normalize to 5m-Median Based on Cohort", outline = F)
points(jitter(rep(1,10), factor = 3),my.data$`5m`, pch = 16, col = c(rep("red",5),rep("blue",5)))
points(jitter(rep(2,10), factor = 3),my.data$`21m`, pch = 16, col = c(rep("red",5),rep("blue",5)))
#points(jitter(rep(3,13), factor = 3),my.data$`12m`, pch = 16, col = c(rep("black",2)))
#points(jitter(rep(4,13), factor = 3),my.data$`25m`, pch = 16, col = c(rep("black",2)))
legend("topright",inset = 0.02,legend=c("cohort 1", "cohort 2"),fill = c("red","blue"))
dev.off()

wilcox.test(my.data$`5m`,my.data$`21m`)

#wilcox.test(my.data$`5m`,my.data$`12m`)

#wilcox.test(my.data$`12m`,my.data$`21m`)

#wilcox.test(my.data$`12m`,my.data$`25m`)

#wilcox.test(my.data$`5m`,my.data$`25m`)

#wilcox.test(my.data$`21m`,my.data$`25m`)

