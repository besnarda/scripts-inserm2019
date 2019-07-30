library(ggplot2)
data <- read.table("WORK/SNP_BX_Ulcerans.txt")
ggplot() + geom_density(aes(data$V1),adjust=1/10)+theme_classic()

data <- read.table("WORK/SNP_CP35_Ulcerans.txt")
ggplot() + geom_density(aes(data$V1),adjust=1/200)+theme_classic()

