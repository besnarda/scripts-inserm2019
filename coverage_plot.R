data <- read.table("/home/t-iris-005/OLD/Mycolactone-project/BAM/S6_coverage.txt")


library(ggplot2)
data2 <- data[data$V3>100,]
ggplot(data=data2, aes(y=log(V3),x=V2)) +geom_point() +  facet_grid(.~V1)

##########
data <- read.table("/home/t-iris-005/OLD/Mycolactone-project/BAM/S7_coverage.txt")
data2 <- data[data$V3>100,]
ggplot(data=data2, aes(y=log(V3),x=V2)) +geom_point() +  facet_grid(.~V1)
                                                                    

##########
data <- read.table("/home/t-iris-005/OLD/Mycolactone-project/BAM/S8_coverage.txt")
data2 <- data[data$V3>100,]
ggplot(data=data2, aes(y=log(V3),x=V2)) +geom_point() +  facet_grid(.~V1)
