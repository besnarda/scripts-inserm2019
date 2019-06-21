
library(ggplot2)
library(ggmap)

#data <- read.csv("0-RAW_DATA/Fichier_bilan_souches.csv",dec=",")
data <- read_xlsx("/home/t-iris-005/0-RAW_DATA/Summary_all_strains_2019.xlsx")
register_google(key = "AIzaSyC0aNyGZX5ZZRzfEI8ybIjz5Do1--lkweE")

summary(data)
data$Latitude <- as.numeric(data$Latitude)
data$Longitude <- as.numeric(data$Longitude)

#subset Benin Nigeria
data <- data[data$Country %in% c("Benin","Nigeria") & data$Cluster != "?",]
ggplot(data=data,aes(x=Div1,fill=Cluster)) + geom_bar()

### plot all strains depending on genotypes
Benin <- c(left = 1.6, right = 4, bottom = 6.2, top = 9)
benin_map <- get_stamenmap(Benin,maptype= "toner-lite") %>% ggmap

benin_map + geom_jitter(aes(y=Latitude,x=Longitude,color=Cluster),data=data)

# sans kouffo
Benin <- c(left = 2.2, right = 2.8, bottom = 6.2, top = 7.5)
benin_map <- get_stamenmap(Benin,maptype= "toner-lite") %>% ggmap

benin_map + geom_jitter(aes(y=Latitude,x=Longitude,color=Cluster),data=data) + facet_wrap(.~Cluster)


benin_map<- get_map(location = c(lon = 2.6, lat = 6.85),
                    color = "color",
                    source = "google",
                    maptype = "terrain",
                    zoom = 9) %>% ggmap()

benin_map + 
  geom_jitter(aes(y=Latitude,x=Longitude,color=Cluster),data=data) +
  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d','#666666',"blue","red"))
benin_map + 
  geom_jitter(aes(y=Latitude,x=Longitude,color=Cluster),data=data) + facet_wrap(.~Cluster)+
  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d','#666666',"blue","red"))


#### stat sur le cameroun

library(readODS)
library(maps)
library(mapdata)

data<- read_ods("0-RAW_DATA/Summary_all_strains_2019.ods")

data$Latitude <- as.numeric(data$Latitude)
cameroon <- data[data$Ref =="Us_cameroun",]
summary(cameroon)

### test de carte avec google map.


map_border <- c(left = 11.8, right = 12.5, bottom = 3.5, top = 4)

map_cameroon <- get_map(location = c(lon = 12.2, lat = 3.8),
                         color = "color",
                         source = "google",
                         maptype = "terrain",
                         zoom = 10) %>% ggmap

map_cameroon + 
  geom_point(data=cameroon,aes(y=Latitude,x=Longitude,color=Cluster))+
  geom_label_repel(data=cameroon,aes(y=Latitude,x=Longitude,color=Cluster,label=`Isolate number`))
