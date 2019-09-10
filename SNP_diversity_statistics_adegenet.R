### calculating statistic for SNP data

## loading libraries
# global
library(tidyverse)
library(readxl)
# specific
library(PopGenome)

snp <- readData("A-GENOME_TO_TREE/2-SNIPPY_CORE/popgenome-vcf",format="VCF")
data <- read_xlsx("/home/t-iris-005/0-RAW_DATA/Summary_all_strains_2019.xlsx")
names <- get.individuals(snp)[[1]]
pop <- left_join(as.data.frame(names),data[data$Dup1_Bad2==0,],by=c("names"="Run")) %>% select(Cluster,names)

# setting population 1
pop1 <- pop[pop$Cluster=="1","names"]
pop2 <- pop[pop$Cluster=="2","names"]
pop3 <- pop[pop$Cluster=="3","names"]
pop4 <- pop[pop$Cluster=="4","names"]
pop5 <- pop[pop$Cluster=="5","names"]
pop6 <- pop[pop$Cluster=="6","names"]
pop7 <- pop[pop$Cluster=="7","names"]
pop8 <- pop[pop$Cluster=="8","names"]
popMuA2 <- pop[pop$Cluster=="Mu_A2","names"]
popKouffo <- pop[pop$Cluster=="Kouffo","names"]

snp <- set.populations(snp,list(pop1,pop2,pop3,pop4,pop5,pop6,pop7,pop8,popMuA2,popKouffo))
snp@populations 
snp <- F_ST.stats(snp)
get.diversity(snp)[,]
get.F_ST(snp)[,]

### setting population to calculate overall diversity
popMuA1 <- pop[pop$Cluster %in% c("1","2","3","4","5","6","7","8","Agogo-1","Kouffo","?"),"names"]
popMuA2 <- pop[pop$Cluster=="Mu_A2","names"]
snp <- set.populations(snp,list(popMuA1,popMuA2))
snp@populations 
snp <- F_ST.stats(snp)
get.diversity(snp)[,]
snp@nuc.diversity.within[1,]
snp@nuc.diversity.between

snp <- neutrality.stats(snp)
get.neutrality(snp)[,]
