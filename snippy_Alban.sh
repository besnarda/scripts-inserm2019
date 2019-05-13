#!/bin/bash

# exporting varaibles and adding emboss to path
export PERL5LIB=/home/t-iris-005/.cpan/build/BioPerl-1.007002-A8KCCK
export PATH=$PATH:/home/t-iris-005/SOFTWARE/vt::/home/t-iris-005/Telechargements/EMBOSS-6.6.0/emboss::/home/t-iris-005/SOFTWARE/minimap2/


##### STEP 1 to launch one by one our analysis

# variables
SAMPLE="/home/t-iris-005/0-RAW_DATA/Genomes_Cameroun/11_27_R1.fastq"
OUTDIR="/home/t-iris-005/A-GENOME_TO_TREE/1-SNIPPY/11_27"
REF="/home/t-iris-005/0-RAW_DATA/References/Agy99_pMUM001_Hard_masked.gbff"
# launch snippy
/home/t-iris-005/SOFTWARE/snippy/bin/snippy --cpus 8 --outdir $OUTDIR --reference $REF --R1 $SAMPLE --R2 ${SAMPLE/_R1/_R2} 

###### STEP 2 to launch snippy for multiple data

# A lancer dans le répertoire on l'on veut créer les fichiers.

for SAMPLE in $(ls /home/t-iris-005/Genomes_Cameroun/*R1.fastq); do
    echo $SAMPLE;
    OUTDIR=$(basename $SAMPLE | sed 's/_R1.fastq//');
    REF="/home/t-iris-005/0-RAW_DATA/References/Agy99_pMUM001_Hard_masked.gbff";   
    # launch snippy
    /home/t-iris-005/snippy/bin/snippy --cpus 8 --outdir $OUTDIR --reference $REF --R1 $SAMPLE --R2 ${SAMPLE/_R1/_R2} ;
done




###### STEP 3 to lauch snippy-core (need previously done snippy analysis)

# !! you need to be in the folder you want the output !!

# variables
SNIPPY_FOLDERS_LIST=$(cat list_383_strains.txt)  # will give all folders in the directory
REF="/home/t-iris-005/0-RAW_DATA/References/Agy99_pMUM001_Hard_masked.gbff"
PREFIX="383_strains"

# launch snippy-core
/home/t-iris-005/SOFTWARE/snippy/bin/snippy-core --prefix $PREFIX  $SNIPPY_FOLDERS_LIST


##### some supplementary stats on snippy-core results
# get SNP name
for tab in *tab; do cat $tab | cut -f 1,2 | sed 's/\t/_/g' | sed '1d' > ${tab}.SNP; done
# compare SNP name between two files (nombre de SNP commun)
comm -12 <(sort 179_strains_no_outgroup.tab.SNP) <(sort 202_strains_no_outgroup.tab.SNP) | wc -l



