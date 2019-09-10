# The objectives of this script is to treat our peak data output
# It will allow me to imporve in dplyr and the use of DNA link function in R with seqinR


library(tidyverse)
library(seqinr)
library(gridExtra)
library(Peptides)

gene_info <- read.table("~/0-RAW_DATA/References/Agy99.gff3",sep="\t",quote = "")
peak_data <- read.table("~/GEXPLORE-MYCO-2019/03b-PEAK_CALLING/myco_bam_peaks.xls",header=TRUE)
burulist <- read.table("~/0-RAW_DATA/Database_Burulist.xls",header=TRUE,sep="\t",quote="")

# some little plots to discover peak data
ggplot(data=peak_data,aes(y=log10(pileup),x=abs_summit))+geom_point()+facet_wrap(.~chr,scales="free_x")
ggplot(data=peak_data,aes(x=length))+geom_histogram(color="black", fill="white")
ggplot(data=peak_data,aes(y=log10(pileup),x=length,color=chr))+geom_point()

# renaming gff3 colnames file 
colnames(gene_info) <- c("sequence","source","feature","start","end","score","strand","phase","attributes")
# keeping only gene and pseudogene data
gene_info <- gene_info %>% filter(feature %in% c("gene","pseudogene"))

# there is some missing entry in database ncbi  compared to database burulist (ex: MUL_0031)

# assign gene or pseudogene to each peak
# we will only do that for peak overlapping with one gene only
# i.e. if a peak overlap two genes we don't write any!

# a function to return gene at a given position
get_gene_at_position <- function(gene_info,position,chr){
  return(gene_info[gene_info$start <= position & gene_info$end >= position & gene_info$sequence ==chr,])
}

# a function to return all genes overlapping with a peak
get_gene_for_peak <- function(gene_info,peak){
  chr <- as.character(peak$chr)
  start <- get_gene_at_position(gene_info,peak$start,chr)
  end <- get_gene_at_position(gene_info,peak$end,chr)
  return(unique(rbind(start,end)))
}

# adding a column indicating the gene in our peak_data
for (i in 1:nrow(peak_data)){
  peak = peak_data[i,]
  tmp_gene <- get_gene_for_peak(gene_info,peak)
  if (nrow(tmp_gene) ==1) {
    tmp_gene <- str_split(string=tmp_gene$attributes,";")[[1]][1]
    tmp_gene <- str_split(tmp_gene,"ID=gene-")[[1]][2]
    peak_data[i,"gene"] <- tmp_gene
  }
}

# Genes with at least 3 peaks
table(as.factor(peak_data$gene))[table(as.factor(peak_data$gene))>=3]
peak_enriched <- merge(peak_data,burulist,by.x="gene",by.y="Accession.number",all.x =TRUE)

# ajouter un graph proportion génome complet selon les catégories. (en multipliant par la taille des gènes)
# Le comparer à notre graph pour voir si la probabilité de tomber 
# dans un gène de tel ou tel type est forte ou pas
palette <- c('#fdb462','#8dd3c7','#bebada','#fb8072','#80b1d3','#ffffb3','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5','#ffed6f')
p1 <- ggplot(data=burulist,aes(x="Expected",fill=BURUlist))+geom_bar()+scale_fill_manual(values = palette)+guides(fill=FALSE)
p2 <- ggplot(data=peak_enriched,aes(fill=BURUlist,x="Observed"))+geom_bar()+scale_fill_manual(values = palette,na.value="darkgrey")
grid.arrange(p1, p2,nrow=1,widths=c(1,2))
ggplot(data=peak_enriched,aes(y=log10(pileup),x=length,color=BURUlist))+geom_point()

##################################################################################################
#### traduction dans les gènes, étude des peptides
Agy99.fa <- read.fasta("~/0-RAW_DATA/References/Agy99.fa")
count=0
for (i in 1:nrow(peak_enriched)){
  
  peak1 <- peak_enriched[i,]
  # on vérifie qu'il y a bien un gène associé
  if (is.na(peak1$gene)){peak1.prot <- NA; next}
  # on récupère la séquence du peak
  peak1.seq <- Agy99.fa[[as.character(peak1$chr)]][peak1$start:peak1$end]
  # on calcule la frame nécessaire pour être comme dans le gène associé
  peak1.gene <- gene_info[grepl(paste0("ID=gene-",peak1$gene),gene_info$attributes),c("start","end","strand")]
  peak1.frame <- (peak1.gene$start - peak1$start) %% 3
  # ça marche
  if ( peak1.gene$strand == "+"){
    peak1.frame <- (peak1.gene$start - peak1$start) %% 3
    peak1.prot <- translate(peak1.seq,frame=peak1.frame)
  } else if (peak1.gene$strand =="-"){
    peak1.frame <- (peak1$end - peak1.gene$end) %% 3
    peak1.prot <- translate(peak1.seq,frame=peak1.frame,sens="R")
  }
  if("*" %in% peak1.prot){ count = count+1;peak1.prot <- NA; next}
  peak_enriched[i,"prot"] <- paste0(peak1.prot,collapse="")
  
  
  
}
count

peak_enriched$hydrophobicity <- hydrophobicity(peak_enriched$prot)



write.table(x=peak_enriched,file="~/GEXPLORE-MYCO-2019/03b-PEAK_CALLING/peak_enriched.xls",quote=FALSE,sep="\t",row.names = FALSE)
# 94 peptide ont un codon stop (pseudogène ou portion hors gene)! reste 1009-94 = 915 peptides restants


########################


#exploration des data avec tidyverse
library(tidyverse)
peak_enriched <- read.table("~/GEXPLORE-MYCO-2019/03b-PEAK_CALLING/peak_enriched.xls",header=TRUE,sep="\t",quote="")

# histogramme de l'hydrophobicité des 913 peptides restants 
# 0.4 c'est suffisament hydrophobe pour être une hélice est 0.7 c'est très hydrophobe (Eisenberg 1984)
peak_enriched %>% filter(hydrophobicity !=0) %>% 
  ggplot(aes(x=hydrophobicity))+ 
  geom_histogram() +
  geom_vline(xintercept = 0.4)+
  geom_vline(xintercept = 0.7)

# hydrophobicité VS pileup avec les catégories en couleurs...
palette <- c('#fdb462','#8dd3c7','#bebada','#fb8072','#80b1d3','#ffffb3','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5','#ffed6f')
peak_enriched %>% filter(hydrophobicity !=0) %>% 
  ggplot(aes(y=log10(pileup),x=hydrophobicity,color=BURUlist))+
  geom_point()+
  scale_color_manual(values=palette)

peak_enriched %>% filter(BURUlist == '"Virulence, detoxification, adaptation"')

peak_enriched %>% ggplot(aes(x=log10(pileup))) + geom_histogram()
