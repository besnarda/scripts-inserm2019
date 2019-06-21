
# package ape is well known for tree visualization
# this script aims at drawing trees from Phyml software
# input is a list of tree in newick format.
# output will be different graph and stat.

### charging libraries
library(ggtree)  # avaible with bioconductor
library(ape)
library(readxl)

####################################################
###### tree with old sequences (184)
####################################################
tree <- read.tree("/home/t-iris-005/A-GENOME_TO_TREE/3-PHYML/184_strains.phylip_phyml_tree_184_strains.phylip.txt")
angers_bilan  <- read.csv("/home/t-iris-005/0-RAW_DATA/Fichier_bilan_souches.csv")

# list of individuals (strains in data)
indiv <- tree$tip.label

# removing outgroup
outgroup_name = c("1172-13","1115-16","856-16","320-7-12","282-5-11","283-8-11","480-14","ERX932498","Reference","270-3-11")
tree_test <- drop.tip(tree,tip=outgroup_name)

# coloring individuals based on previously identified genotypes
Group_geno = list(as.character(angers_bilan[angers_bilan$Genotype==1,"MU.strain.rename.angers"]),
                  as.character(angers_bilan[angers_bilan$Genotype==2,"MU.strain.rename.angers"]),
                  as.character(angers_bilan[angers_bilan$Genotype==3,"MU.strain.rename.angers"]),
                  as.character(angers_bilan[angers_bilan$Genotype==4,"MU.strain.rename.angers"]),
                  as.character(angers_bilan[angers_bilan$Genotype==5,"MU.strain.rename.angers"]),
                  as.character(angers_bilan[angers_bilan$Genotype==6,"MU.strain.rename.angers"]),
                  as.character(angers_bilan[angers_bilan$Genotype==7,"MU.strain.rename.angers"]),
                  as.character(angers_bilan[angers_bilan$Genotype==8,"MU.strain.rename.angers"]))
tree_test <- groupOTU(tree_test,Group_geno)

# plotting the tree
ggtree(tree_test,aes(color=group),layout="circular") +theme_tree2() + geom_tiplab2(size=2.5) 


####################################################
###### tree with old sequences + new from 2014-2016
####################################################

# getting data
tree <- read.tree("/home/t-iris-005/OLD/PhyloTrees/Trees_Aout2018/final_trees/Aout2018_750coll.aln_phyml_tree_final_b1000_core.txt")
tree <- read.tree("/home/t-iris-005/A-GENOME_TO_TREE/3-PHYML/208_strains.phylip_phyml_tree_208_strains.phylip.txt")
angers_bilan  <- read.csv("/home/t-iris-005/0-RAW_DATA/Fichier_bilan_souches.csv")
country_bilan <- read.csv("/home/t-iris-005/0-RAW_DATA/strains_country.txt")

# list of individuals (strains in data)
indiv <- tree$tip.label

# removing outgroup
outgroup_name = c("1172-13","1115-16","856-16","320-7-12","282-5-11","283-8-11","480-14","ERX932498","Reference","270-3-11")
tree_test <- drop.tip(tree,tip=outgroup_name)

# coloring individuals based on previously identified genotypes
Group_geno = list(as.character(angers_bilan[angers_bilan$Genotype==1,"MU.strain.rename.angers"]),
              as.character(angers_bilan[angers_bilan$Genotype==2,"MU.strain.rename.angers"]),
              as.character(angers_bilan[angers_bilan$Genotype==3,"MU.strain.rename.angers"]),
              as.character(angers_bilan[angers_bilan$Genotype==4,"MU.strain.rename.angers"]),
              as.character(angers_bilan[angers_bilan$Genotype==5,"MU.strain.rename.angers"]),
              as.character(angers_bilan[angers_bilan$Genotype==6,"MU.strain.rename.angers"]),
              as.character(angers_bilan[angers_bilan$Genotype==7,"MU.strain.rename.angers"]),
              as.character(angers_bilan[angers_bilan$Genotype==8,"MU.strain.rename.angers"]))
tree_test <- groupOTU(tree_test,Group_geno)

# plotting the tree
ggtree(tree_test,aes(color=group),layout="circular") +theme_tree2() + geom_tiplab2(size=2.5) 
### fig with country 
Group_geno <- list(as.character(country_bilan[country_bilan$Country=="Nigeria","Strain"]),
                   as.character(country_bilan[country_bilan$Country=="Benin","Strain"]))
tree_test <- groupOTU(tree_test,Group_geno)
# plotting the tree
ggtree(tree_test,aes(color=group),layout="circular") +theme_tree2() + geom_tiplab2(size=2.5) +
  scale_color_manual(values=c('#7570b3','#a6761d'),
                     name= "Legend", guide = "legend", labels= c("Nigeria","Benin"))+
  theme(legend.position = c(0.2,0.2), 
        legend.title = element_blank(), # no title
        legend.key = element_blank()) # no keys


####################################################
###### tree with all our sequences + data from vanhoote 2019 + data cameroun + data ghana
####################################################

# getting data
tree <- read.tree("/home/t-iris-005/A-GENOME_TO_TREE/3-PHYML/425_strains.phylip_phyml_tree_425_strains.phylip.txt")
angers_bilan  <- read.csv("/home/t-iris-005/0-RAW_DATA/Fichier_bilan_souches.csv")
vander_bilan  <- read.csv("/home/t-iris-005/OLD/snippy_WGHM_Feb2019/ID_table.csv")
country_bilan <- read.csv("/home/t-iris-005/0-RAW_DATA/strains_country.txt")

# list of individuals (strains in data)
indiv <- tree$tip.label
tree$tip.label <- gsub("^n","",x= tree$tip.label)

# removing outgroup
outgroup_name = c("10-381","1172-13","1115-16","856-16","320-7-12","282-5-11","283-8-11","480-14","ERX932498","270-3-11")
tree_test <- drop.tip(tree,tip=outgroup_name)
# tree_test <- tree

##### fig2. coloring individuals based on previously identified genotypes (ANGERS)
Group_geno = list(as.character(angers_bilan[angers_bilan$Genotype==1,"MU.strain.rename.angers"]),
                  as.character(angers_bilan[angers_bilan$Genotype==2,"MU.strain.rename.angers"]),
                  as.character(angers_bilan[angers_bilan$Genotype==3,"MU.strain.rename.angers"]),
                  as.character(angers_bilan[angers_bilan$Genotype==4,"MU.strain.rename.angers"]),
                  as.character(angers_bilan[angers_bilan$Genotype==5,"MU.strain.rename.angers"]),
                  as.character(angers_bilan[angers_bilan$Genotype==6,"MU.strain.rename.angers"]),
                  as.character(angers_bilan[angers_bilan$Genotype==7,"MU.strain.rename.angers"]),
                  as.character(angers_bilan[angers_bilan$Genotype==8,"MU.strain.rename.angers"]))
tree_test <- groupOTU(tree_test,Group_geno)

# plotting the tree
ggtree(tree_test,aes(color=group)) +theme_tree2() + geom_tiplab(size=2.5) 

##### fig3. coloring individuals based on vandernoote information

Group_geno <- list(as.character(vander_bilan[vander_bilan$BAPS_group=="BAPS_1","Experiment"]),
                   as.character(vander_bilan[vander_bilan$BAPS_group=="BAPS_2","Experiment"]),
                   as.character(vander_bilan[vander_bilan$BAPS_group=="BAPS_3","Experiment"]),
                   as.character(vander_bilan[vander_bilan$BAPS_group=="BAPS_4","Experiment"]),
                   as.character(vander_bilan[vander_bilan$BAPS_group=="BAPS_5","Experiment"]),
                   as.character(vander_bilan[vander_bilan$BAPS_group=="BAPS_6","Experiment"]))
tree_test <- groupOTU(tree_test,Group_geno)
# plotting the tree
ggtree(tree_test,aes(color=group)) +theme_tree2() + geom_tiplab(size=2.5)

#### fig4. coloring depending on the country
Group_geno <- list(as.character(country_bilan[country_bilan$Country=="DRC","Strain"]),
                   as.character(country_bilan[country_bilan$Country=="Nigeria","Strain"]),
                   as.character(country_bilan[country_bilan$Country=="Cameroun","Strain"]),
                   as.character(country_bilan[country_bilan$Country=="Ghana","Strain"]),
                   as.character(country_bilan[country_bilan$Country=="Angola","Strain"]),
                   as.character(country_bilan[country_bilan$Country=="Benin","Strain"]),
                   as.character(country_bilan[country_bilan$Country=="Congo","Strain"]))
tree_test <- groupOTU(tree_test,Group_geno)
# plotting the tree
ggtree(tree_test,aes(color=group)) + geom_tiplab(size=2.5) +
  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d','#666666'),
                     name= "Legend", guide = "legend", labels= c("No Data","RDC","Nigeria","Cameroun","Ghana","Angola","Benin","Congo"))+
  theme(legend.position = c(0.5,0.5), 
        legend.title = element_blank(), # no title
        legend.key = element_blank()) # no keys


pdf(file="/home/t-iris-005/A-GENOME_TO_TREE/4-R_PLOT/425_strains_country.pdf",width = 100,height =200)
ggtree(tree_test,aes(color=group))  + geom_tiplab(size=10)+
  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d','#666666'),
                     name= "Legend", guide = "legend", labels= c("No Data","RDC","Nigeria","Cameroun","Ghana","Angola","Benin","Congo"))+
  theme(legend.position = c(0.5,0.5),legend.text=element_text(size=100),
        legend.title = element_blank(), # no title
        legend.key = element_blank()) 

dev.off()


####################################################
###### tree for Cameroun strains
####################################################
tree <- read.tree("/home/t-iris-005/A-GENOME_TO_TREE/3-PHYML/22_strains.phylip_phyml_tree_22_strains.txt")

outgroup_name = c("320-7-12","Reference")
tree <- root(phy = tree,outgroup_name,resolve.root=T)
tree_test <- drop.tip(tree,tip=outgroup_name)

ggtree(tree_test) +theme_tree2() + geom_tiplab(size=5)

####################################################
###### tree for Benin and nigeria strains
####################################################

tree <- read.tree("/home/t-iris-005/A-GENOME_TO_TREE/3-PHYML/benin_nigeria_strains.phylip_phyml_tree_benin_nigeria_strains.txt")
info <- read_xlsx("/home/t-iris-005/0-RAW_DATA/Summary_all_strains_2019.xlsx")

# on crée un data.frame avec toutes les infos que l'on connait sur les points
info <- data.frame(info)
info <- merge(data.frame(tree$tip.label),info,by.x="tree.tip.label",by.y="Run")

# rooting and deleting outgroup Here we use Mu_A2!
outgroup = as.character(info[is.na(info$Lineage)==FALSE & info$Lineage == "Mu_A2","tree.tip.label"])

# pour décider ce qui fait aprtie de l'outgroup!
# tree_test <- drop.tip(tree,tip=outgroup)
# ggtree(tree) +theme_tree2() + geom_tiplab(size=2.5)

# je n'arrive pas à faire le graph avec l'ensemble des Mu_A2, on ne prend qu'un pour faire l'outgroup
tree <- root(tree,outgroup[1],resolve.root=T)
tree_test <- drop.tip(tree,tip=outgroup)

# plot with colour depending of Country
pdf(file="/home/t-iris-005/A-GENOME_TO_TREE/4-R_PLOT/Benin_Nigeria_country.pdf")

ggtree(tree_test) %<+% info +
  theme_tree2() +
  geom_tiplab(size=0.8,aes(label=Isolate.number)) +
  geom_tippoint(size =0.5,aes(color=Country))+
  theme(legend.position="right")

dev.off()

# plot with colour depending of Cluster
pdf(file="/home/t-iris-005/A-GENOME_TO_TREE/4-R_PLOT/benin_nigeria_8geno.pdf")

ggtree(tree_test,aes(color=Cluster)) %<+% info +
  theme_tree2() +
  geom_tiplab(size=0.8,aes(label=Isolate.number)) +
  theme(legend.position="right")

dev.off()
####################################################
###### tree for Cameroun strains
###################################################

tree <- read.tree("/home/t-iris-005/A-GENOME_TO_TREE/3-PHYML/cameroun_strains.phylip_phyml_tree_cameroun_strains.txt")

# on crée un data.frame avec toutes les infos que l'on connait sur les points
info <- read_xlsx("/home/t-iris-005/0-RAW_DATA/Summary_all_strains_2019.xlsx")
info <- data.frame(info)
info <- merge(data.frame(tree$tip.label),info,by.x="tree.tip.label",by.y="Run")

# rooting and deleting outgroup Here we use Mu_A2!
outgroup = as.character(info[is.na(info$Lineage)==FALSE & info$Lineage == "Mu_A2","tree.tip.label"])
tree <- root(tree,outgroup,resolve.root=T)
tree_test <- drop.tip(tree,tip=outgroup)

# plot with colour depending of village. (lot of village so not really easy...)
ggtree(tree_test) %<+% info +
  theme_tree2() +
  geom_tiplab(size=2,aes(label=Isolate.number)) +
  geom_tippoint(aes(color=Div3))+
  theme(legend.position="right")


  
