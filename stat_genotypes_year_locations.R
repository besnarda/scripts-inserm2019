
library(ggplot2)
library(ggmap)
library(readxl)
library(ggrepel)
library(ggforce)

#data <- read.csv("0-RAW_DATA/Fichier_bilan_souches.csv",dec=",")
data <- read_xlsx("/home/t-iris-005/0-RAW_DATA/Summary_all_strains_2019.xlsx")
register_google(key = "AIzaSyC0aNyGZX5ZZRzfEI8ybIjz5Do1--lkweE")

summary(data)
data$Latitude <- as.numeric(data$Latitude)
data$Longitude <- as.numeric(data$Longitude)

my_color_scale = c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d','#666666',"darkgreen","blue")

#subset Benin Nigeria
data <- data[data$Country %in% c("Benin","Nigeria") & data$Cluster != "?",]
ggplot(data=data,aes(x=Div1,fill=Cluster)) + geom_bar()

### plot all strains depending on genotypes
Benin <- c(left = 1.6, right = 4, bottom = 6.2, top = 9)
benin_map <- get_stamenmap(Benin,maptype= "toner-lite") %>% ggmap

benin_map + geom_jitter(aes(y=Latitude,x=Longitude,color=Cluster),data=data) +scale_color_manual(values=my_color_scale)

# sans kouffo
Benin <- c(left = 2.2, right = 2.8, bottom = 6.2, top = 7.5)
benin_map <- get_stamenmap(Benin,maptype= "toner-lite") %>% ggmap

benin_map + geom_point(aes(y=Latitude,x=Longitude,color=Cluster),size = 2,data=data[data$Cluster %in% c("Kouffo","Mu_A2")==FALSE,]) +scale_color_manual(values=my_color_scale)
benin_map + geom_point(aes(y=Latitude,x=Longitude,color=Cluster),size = 2,position = position_jitter(w = 0.03, h = 0.03),data=data[data$Cluster %in% c("Kouffo","Mu_A2")==FALSE,]) +scale_color_manual(values=my_color_scale)


benin_map<- get_map(location = c(lon = 2.6, lat = 6.85),
                    color = "color",
                    source = "google",
                    maptype = "terrain",
                    zoom = 9) %>% ggmap()

benin_map +
  geom_jitter(aes(y=Latitude,x=Longitude,color=Ref), data=data)
benin_map +
  geom_jitter(aes(y=Latitude,x=Longitude,color=Cluster),data=data) +
  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d','#666666',"blue","red"))
benin_map +
  geom_jitter(aes(y=Latitude,x=Longitude,color=Cluster),data=data) + facet_wrap(.~Cluster)+
  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d','#666666',"blue","red"))

# plot avec juste un génogroupe et des labels
benin_map + ggrepel::geom_label_repel(aes(y=Latitude,x=Longitude,color=Cluster,label=`Isolate number`),data=data[data$Cluster=="8",])

#### ajout des clusters détectés par satscan!

circles <- data.frame(lon=c(6.95,6.58,6.93,7.06),lat=c(3.51,2.55,2.27,2.62),r=c(85.36/100,10.55/100,24.83/100,18.61/100),label=c("Nigeria","Sud_oueme","Nord_oueme","Plateau"))
benin_map + geom_circle(data=circles,aes(y0=lon,x0=lat,r=r)) +
  geom_jitter(aes(y=Latitude,x=Longitude,color=Cluster),data=data)

#########################################################
#### map sur le cameroun
#########################################################

cameroon <- data[data$Country =="Cameroon" & data$Cluster!="BAPS-2" & data$Dup1_Bad2==0,]
summary(cameroon)

### test de carte avec google map.


map_border <- c(left = 11.8, right = 12.5, bottom = 3.5, top = 4)

map_cameroon <- get_map(location = c(lon = 12.2, lat = 3.8),
                         color = "color",
                         source = "google",
                         maptype = "terrain",
                         zoom = 9) %>% ggmap

map_cameroon +
  geom_point(data=cameroon,aes(y=Latitude,x=Longitude,color=as.factor(Cluster)))+
  geom_label_repel(data=cameroon,aes(y=Latitude,x=Longitude,color=Cluster,label=`Isolate number`))

map_cameroon +
  geom_jitter(data=cameroon,aes(y=Latitude,x=Longitude,color=as.factor(Cluster)))
