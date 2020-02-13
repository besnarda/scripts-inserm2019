
#--------------------------------------------------------------------------------------------#
# 		SCRIPT : PLOT A PHYML TREE
#					With a phyml tree and supplementary information. Draw the tree
#
#--------------------------------------------------------------------------------------------#

#### Arguments management with docopt
###################################################################

'plot_phyml_tree.
this is some additional text
Usage:
  plot_phyml_tree.R <tree_file> <info_file> <output_file> [--remove=<remove_file>] [--outgroup=<outgroup_file>] [--color=<column_name>]
  plot_phyml_tree.R (-h | --help)
  plot_phyml_tree.R --version

Options:
  -h --help     Show this screen.
  --version     Show  version.
  --remove=<remove_file>  a file with a list of strains you want to remove from the final tree. [default: NA]
  --outgroup=<outgroup_file>  a file with the list of outgroup strains. [default: Reference]
  --color=<column_name>   the name of the column to use as factor. [default: Cluster]
  

' -> doc

library(docopt)
version = 'plot_phyml_tree 0.1'
print(paste0("This is ",version))
arguments <- docopt(doc, version = version)

# test with 
# arguments <- docopt(doc,args = c("phyml/425_strains_phyml_tree.txt","data/summary_strains.xlsx","--color Country"),version = "plot_phyml_tree 0.1")
# print(arguments)

##################################################

tree_file = arguments$tree_file
info_file = arguments$info_file
outgroup_file = arguments$outgroup
remove_file = arguments$remove
column_name = arguments$color
output_file = arguments$output_file


#### loading librairies

suppressWarnings(suppressMessages(library(ggtree)))  # avaible with bioconductor
suppressWarnings(suppressMessages(library(ape)))
library(readxl)
library(ggplot2)

#### importing data and merging
tree <- read.tree(tree_file)
info <- read_xlsx(info_file)

info <- data.frame(info)
info <- merge(data.frame(tree$tip.label),info,by.x="tree.tip.label",by.y="Run")

#### defining the outgroup
if (outgroup_file != "Reference"){
  outgroup_list = read.csv(outgroup_file)
} else {
  outgroup_list = "Reference"
}



#### defining the remove group
if (remove_file == "NA"){
  remove=""
} else {
  remove=read.csv(remove_file)
}

tree <- drop.tip(tree,tip=remove)

#### defining the outgroup
outgroup = as.character(info[
  (is.na(info$Isolate.number)==FALSE & info$Isolate.number %in% outgroup_list) |
    (is.na(info$tree.tip.label)==FALSE & info$tree.tip.label %in% outgroup_list) |
    (is.na(info$Cluster)==FALSE & info$Cluster %in% outgroup_list),"tree.tip.label"])
tree <- root(tree,outgroup,resolve.root=T)
# drop outgroup
tree_test <- drop.tip(tree,tip=outgroup)

my_color_scale = c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d','#666666',"darkgreen","blue")

# getting the best height for svg file
fitted_tiplab_size = 160/length(tree_test$tip.label)
fitted_tippoint_size = 66/length(tree_test$tip.label)

p <- ggtree(tree_test) %<+% info +
  theme_tree2() +
  geom_tiplab(size=fitted_tiplab_size) +
  geom_tippoint(aes(color=get(column_name)),size=fitted_tippoint_size)+
  theme(legend.position="right")+
  scale_color_manual(values=my_color_scale,na.value="red")+
  labs(color = column_name)+
  guides(colour = guide_legend(override.aes = list(size=2)))

# correction of a visual bug that didn't write the whole name of our samples
# so we expand the x-axis limit
p <- p + xlim(NA, max(p$data$x)*1.2)

ggsave(file=output_file, plot=p , width = 5, height = 7)

