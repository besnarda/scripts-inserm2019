### use of SNP distance to compare with position (draw lines between exact pair)


library(ggplot2)
library(ggmap)
library(readxl)
library(ggforce)
library(vcfR)
library(adegenet)
library(ape)
library(reshape2)
library(dplyr)

data <- read_xlsx("/home/t-iris-005/0-RAW_DATA/Summary_all_strains_2019.xlsx")
data$Latitude <- as.numeric(data$Latitude)
data$Longitude <- as.numeric(data$Longitude)

###############################################################
###### plot for benin nigeria map 

# subseting data
data <- data[data$Country %in% c("Benin","Nigeria") & data$Cluster != "?",]

# read SNP information
SNP <- read.vcfR("A-GENOME_TO_TREE/2-SNIPPY_CORE/benin_nigeria_strains.vcf")
# SNP <- read.vcfR("A-GENOME_TO_TREE/2-SNIPPY_CORE/425_strains.vcf")
geno <- vcfR2genind(SNP)
# calculating pairwise distance
distmatrix <- as.matrix(dist.gene(geno$tab))
distdf <- melt(distmatrix, varnames = c("strain1", "strain2"))
group0 <- distdf %>% filter(value <= 2 ) %>% filter(strain1 != strain2)
# getting position and cluster of points
group0 <- merge(group0,data[,c("Run","Latitude","Longitude","Cluster")],by.x="strain1",by.y="Run")
group0 <- merge(group0,data[,c("Run","Latitude","Longitude","Cluster")],by.x="strain2",by.y="Run")
# getting closest point
group_close <- distdf %>% filter(value < 100) %>% filter(strain1 != strain2) %>% arrange(value) %>% distinct(strain1, .keep_all=TRUE)
group_close <- merge(group_close,data[,c("Run","Latitude","Longitude","Cluster")],by.x="strain1",by.y="Run")
group_close <- merge(group_close,data[,c("Run","Latitude","Longitude","Cluster")],by.x="strain2",by.y="Run")

# charge map
register_google(key = "AIzaSyC0aNyGZX5ZZRzfEI8ybIjz5Do1--lkweE")
benin_map<- get_map(location = c(lon = 2.6, lat = 6.85),
                    color = "color",
                    source = "google",
                    maptype = "terrain",
                    zoom = 8) %>% ggmap()

# here to change colors!!
my_color_scale = c("1"='#1b9e77',"2"='#d95f02',"3"='#7570b3',"4"='#e7298a',"5"='#66a61e',"6"='#e6ab02',"7"='#a6761d',"8"='#666666',"Kouffo"="blue","Mu_A2"="red","?"="black")
# actually plotting the full map 'lines + point'
benin_map +
  geom_segment(data=group0,aes(y=Latitude.x ,x=Longitude.x ,yend=Latitude.y  ,xend=Longitude.y,color=Cluster.x))+
  geom_point(aes(y=Latitude,x=Longitude,color=Cluster),data=data) +
  scale_color_manual(values=my_color_scale)
# one map per cluster
clust = "?"
benin_map +
  geom_segment(data=group0[group0$Cluster.x==clust,],aes(y=Latitude.x ,x=Longitude.x ,yend=Latitude.y  ,xend=Longitude.y))+
  geom_point(aes(y=Latitude,x=Longitude,color=as.numeric(`Isolation year`)),data=data[data$Cluster==clust,])+
  scale_colour_gradient(low="yellow",high="red",limit=c(1995,2017))+
  labs(color='Isolation Year') + xlab("Longitude") + ylab("Latitude") +
  ggtitle(paste("Genetic proximity (<=1SNP) between strains of cluster",clust))
# plot only with lines
benin_map +
  geom_segment(data=group0,aes(y=Latitude.x ,x=Longitude.x ,yend=Latitude.y  ,xend=Longitude.y,color=Cluster.x))+
  scale_color_manual(values=my_color_scale)
# plot closest point
benin_map +
  geom_segment(data=group_close,aes(y=Latitude.x ,x=Longitude.x ,yend=Latitude.y  ,xend=Longitude.y,color=Cluster.x))+
  scale_color_manual(values=my_color_scale)+facet_wrap(.~Cluster.x)


