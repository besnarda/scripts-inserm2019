
# package ape is well known for tree visualization
# this script aims at drawing trees from Phyml software
# input is a list of tree in newick format.
# output will be different graph and stat.

### charging libraries
library(ggtree)  # avaible with bioconductor

####################################################
###### tree with old sequences (184)
####################################################
tree <- read.tree("/home/t-iris-005/A-GENOME_TO_TREE/3-PHYML/184_strains.phylip_phyml_tree_184_strains.phylip.txt")
angers_bilan  <- read.csv("/home/t-iris-005/OLD/QGIS_clusters/MAP_clusters/Fichier_bilan_souches.csv")

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
ggtree(tree_test,aes(color=group)) +theme_tree2() + geom_tiplab(size=2.5) 


####################################################
###### tree with old sequences + new from 2014-2016
####################################################

# getting data
tree <- read.tree("/home/t-iris-005/OLD/PhyloTrees/Trees_Aout2018/final_trees/Aout2018_750coll.aln_phyml_tree_final_b1000_core.txt")
tree <- read.tree("/home/t-iris-005/A-GENOME_TO_TREE/3-PHYML/208_strains.phylip_phyml_tree_208_strains.phylip.txt")
angers_bilan  <- read.csv("/home/t-iris-005/OLD/QGIS_clusters/MAP_clusters/Fichier_bilan_souches.csv")

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
ggtree(tree_test,aes(color=group)) +theme_tree2() + geom_tiplab(size=2.5) 

####################################################
###### tree with all our sequences + data from vanhoote 2019
####################################################

# getting data
tree <- read.tree("/home/t-iris-005/OLD/PhyloTrees/Trees_Feb2019/core.aln.phylip.Feb2019_phyml_tree_Feb2019_100.txt")
angers_bilan  <- read.csv("/home/t-iris-005/OLD/QGIS_clusters/MAP_clusters/Fichier_bilan_souches.csv")
vander_bilan  <- read.csv("/home/t-iris-005/OLD/snippy_WGHM_Feb2019/ID_table.csv")

# list of individuals (strains in data)
indiv <- tree$tip.label
tree$tip.label <- gsub("^n","",x= tree$tip.label)

# removing outgroup
outgroup_name = c("1172-13","1115-16","856-16","320-7-12","282-5-11","283-8-11","480-14","ERX932498","Reference","270-3-11")
tree_test <- drop.tip(tree,tip=outgroup_name)

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


####################################################
###### tree with all our sequences + data from vanhoote 2019 + data cameroun + data ghana
####################################################

# getting data
tree <- read.tree("/home/t-iris-005/A-GENOME_TO_TREE/3-PHYML/401_strains.phylip_phyml_tree.txt")
tree <- read.tree("/home/t-iris-005/A-GENOME_TO_TREE/3-PHYML/425_strains.phylip_phyml_tree_425_strains.phylip.txt")
angers_bilan  <- read.csv("/home/t-iris-005/OLD/QGIS_clusters/MAP_clusters/Fichier_bilan_souches.csv")
vander_bilan  <- read.csv("/home/t-iris-005/OLD/snippy_WGHM_Feb2019/ID_table.csv")
country_bilan <- read.csv("/home/t-iris-005/0-RAW_DATA/strains_country.txt")

# list of individuals (strains in data)
indiv <- tree$tip.label
tree$tip.label <- gsub("^n","",x= tree$tip.label)

# removing outgroup
outgroup_name = c("10-381","1172-13","1115-16","856-16","320-7-12","282-5-11","283-8-11","480-14","ERX932498","Reference","270-3-11")
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


pdf(file="/home/t-iris-005/A-GENOME_TO_TREE/4-R_PLOT/401_strains_country.pdf",width = 100,height =200)
ggtree(tree_test,aes(color=group))  + geom_tiplab(size=10)+
  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d','#666666'),
                     name= "Legend", guide = "legend", labels= c("No Data","RDC","Nigeria","Cameroun","Ghana","Angola","Benin","Congo"))+
  theme(legend.position = c(0.5,0.5),legend.text=element_text(size=100),
        legend.title = element_blank(), # no title
        legend.key = element_blank()) 

dev.off()
