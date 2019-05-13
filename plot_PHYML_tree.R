
# package ape is well known for tree visualization
# this script aims at drawing trees from Phyml software
# input is a list of tree in newick format.
# output will be different graph and stat.

### charging libraries
library(ggtree)  # avaible with bioconductor

####################################################
###### tree with old sequences + new from 2014-2016
####################################################

# getting data
tree <- read.tree("/home/t-iris-005/PhyloTrees/Trees_Aout2018/final_trees/Aout2018_750coll.aln_phyml_tree_final_b1000_core.txt")
angers_bilan  <- read.csv("/home/t-iris-005/QGIS_clusters/MAP_clusters/Fichier_bilan_souches.csv")

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
tree <- read.tree("/home/t-iris-005/PhyloTrees/Trees_Feb2019/core.aln.phylip.Feb2019_phyml_tree_Feb2019_100.txt")
angers_bilan  <- read.csv("/home/t-iris-005/QGIS_clusters/MAP_clusters/Fichier_bilan_souches.csv")
vander_bilan  <- read.csv("/home/t-iris-005/snippy_WGHM_Feb2019/ID_table.csv")

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
