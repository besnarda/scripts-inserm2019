
library(ggplot2)
library(dplyr)

database <- read.csv("~/0-RAW_DATA/Database_Burulist.xls",sep="\t")
control <- read.csv("GEXPLORE-MYCO-2019/03-PEAK_CALLING/control_peak_annotated.txt",sep="\t")
myco <- read.csv("GEXPLORE-MYCO-2019/03-PEAK_CALLING/myco_peak_annotated.txt",sep="\t")


control <- merge(control, database, by.x="gene_name",by.y="Accession.number")
myco <- merge(myco, database, by.x="gene_name",by.y="Accession.number")

# plotting all gene
myco %>% filter(pileup >= 500) %>%  ggplot(aes(x=BURUlist))+geom_bar()
control %>% filter(pileup >= 500) %>%  ggplot(aes(x=BURUlist))+geom_bar()

myco %>% filter(pileup >= 1000) %>%  ggplot(aes(x=BURUlist))+
  geom_bar() +
  coord_flip()
control %>% filter(pileup >= 1000) %>%  ggplot(aes(x=BURUlist))+
  geom_bar() +
  coord_flip() 

myco %>% ggplot(aes(length))+geom_histogram()+xlim(0,2500)
control %>% ggplot(aes(length))+geom_histogram()+xlim(0,2500)
both <- myco[myco$gene_name %in% control$gene_name,]
myco_alone <- myco[myco$gene_name %in% control$gene_name ==FALSE,]

both %>% filter(pileup >= 1000) %>%  ggplot(aes(x=BURUlist))+
  geom_bar() +
  coord_flip()

myco_alone  %>%  ggplot(aes(x=BURUlist))+
  geom_bar() +
  coord_flip()

myco_alone %>% filter(grepl("carrier",gene_function))
