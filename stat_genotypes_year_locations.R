
library(ggplot2)
library(ggmap)

data <- read.csv("0-RAW_DATA/Fichier_bilan_souches.csv",dec=",")

summary(data)

ggplot(data=data,aes(x=department,fill=Genotype)) + geom_bar()
ggplot(data=data,aes(x=Genotype,fill=department)) + geom_bar()
ggplot(data=data,aes(x=as.factor(year),fill=department)) + geom_bar()

ggplot(data=data,aes(y=X,x=Y,color=Genotype)) + geom_jitter()


### plot all strains depending on genotypes

Benin <- c(left = 2.3, right = 4.3, bottom = 6.4, top = 7.3)
benin_map <- get_stamenmap(Benin,maptype= "toner-lite") %>% ggmap

benin_map + geom_jitter(aes(y=X,x=Y,color=Genotype),data=data) + facet_wrap(.~Genotype)

### plot only benin strains
Benin <- c(left = 2.3, right = 2.9, bottom = 6.4, top = 7.3)
benin_map <- get_stamenmap(Benin,maptype= "toner-lite") %>% ggmap

data_benin_only <- data[is.na(data$department) == FALSE & data$department!="nigeria",]

benin_map + geom_jitter(aes(y=X,x=Y,color=Genotype),data=data_benin_only) + facet_wrap(.~Genotype)


