################ Alban BESNARD ########### albanbesnard@hotmail.fr

#########################################################################################
######## Script shell pour faire des arbres phylogénétiques          ####################
######## à partir de nos données et celles de la littérature         ####################
######## GENOMES_TO_TREE_2019                                         ####################
#########################################################################################

#########################################################################################
###                                   SUMMARY                                         ###
#########################################################################################
###################
###################  00. Raw DATA
###################     0.1 getting data with sra toolkit
###################     0.2 checking data with fastqc
###################     0.3 cleaning adaptators if necessary
###################
###################  01. calling SNP with Snippy
###################     1.1 calling 1 by 1
###################     1.2 calling every sample 
###################
###################  02. Get core SNP with Snippy
###################
###################  03. Make trees with phyml (on cluster)
###################
###################  04. Plotting with R
###################
###################  05. some test with BEAST (time component)
###################
##########################################################################################


#-----------------------------------------------------------------------------------------
#################### STEP 00: Raw data ########################
#-----------------------------------------------------------------------------------------

# On va dans le working directory
cd /home/t-iris-005/A-GENOMES_TO_TREE

# on va dans le dossier raw_data 
cd 0-RAW_DATA

#### 0.1 getting data with sra
############################################################

# make a list of all needed file to fetch (column Run)
# this takes a long time!!
# be aware that fasterq create a cache in your ncbi folder!
# Ask for my Summary_all_strains_2019 file containing multiple information
# aggregated from publications over the time (mostly vandelannoote)


# first one is for benin and nigeria from vandelannote et al. 2017 and others
for SRR in $(cat Run_to_fetch1.txt); do
echo $SRR;
~/SOFTWARE/sratoolkit.2.9.6-1-ubuntu64/bin/fasterq-dump $SRR;
done

# second one is for Cameroun from vandelanoote et al. 2017
for SRR in $(cat Run_to_fetch2.txt); do
echo $SRR;
~/SOFTWARE/sratoolkit.2.9.6-1-ubuntu64/bin/fasterq-dump $SRR;
done

# third one is all data from vandelanoote et al. 2019
for SRR in $(cat Run_to_fetch3.txt); do
echo $SRR;
~/SOFTWARE/sratoolkit.2.9.6-1-ubuntu64/bin/fasterq-dump $SRR;
done

# forth is cameroun data from Bolz et al. 2015 (82 acc)
for SRR in $(cat Run_to_fetch4.txt); do
echo $SRR;
~/SOFTWARE/sratoolkit.2.9.6-1-ubuntu64/bin/fasterq-dump $SRR;
done

# fifth is cameroun data from vandelanoote et al. 2017 (duplicates!)
for SRR in $(cat Run_to_fetch5.txt); do
echo $SRR;
~/SOFTWARE/sratoolkit.2.9.6-1-ubuntu64/bin/fasterq-dump $SRR;
done

# we can zip data to gain space!
for fastq in $(ls *.fastq); do
echo $fastq;
gzip $fastq;
done 

#### 0.2 On vérifie la qualité des données avec un FASTQC
##########################################################

mkdir FASTQC
for i in *1.fastq.gz; do
/home/t-iris-005/SOFTWARE/FastQC/fastqc $i -o FASTQC
done

## Seems to be ok, No need to clean ?


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
 ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
 ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

#-----------------------------------------------------------------------------------------
#################### STEP 01: Calling SNP                         ########################
#-----------------------------------------------------------------------------------------


# exporting varaibles and adding emboss to path
export PERL5LIB=/home/t-iris-005/.cpan/build/BioPerl-1.007002-A8KCCK
export PATH=$PATH:/home/t-iris-005/SOFTWARE/vt::/home/t-iris-005/Telechargements/EMBOSS-6.6.0/emboss::/home/t-iris-005/SOFTWARE/minimap2/


##### 1.1 to launch one by one our analysis

# variables
SAMPLE="/home/t-iris-005/A-GENOME_TO_TREE/0-RAW_DATA/SRR3224040_1.fastq.gz"
OUTDIR="/home/t-iris-005/A-GENOME_TO_TREE/1-SNIPPY/SRR3224040"
REF="/home/t-iris-005/0-RAW_DATA/References/Agy99_pMUM001_Hard_masked.gbff"
# launch snippy
/home/t-iris-005/SOFTWARE/snippy/bin/snippy --cpus 8 --outdir $OUTDIR --reference $REF --R1 $SAMPLE --R2 ${SAMPLE/_R1/_R2}  --mincov 3

###### 1.2 to launch snippy for multiple data

# A lancer dans le répertoire on l'on veut créer les fichiers.
FILE_LIST=$(ls /home/t-iris-005/A-GENOME_TO_TREE/0-RAW_DATA/SR*.fastq.gz)
for SAMPLE in $FILE_LIST; do
    echo $SAMPLE;
    OUTDIR=$(basename $SAMPLE | sed 's/_R1.fastq.gz//');
    REF="/home/t-iris-005/0-RAW_DATA/References/Agy99_pMUM001_Hard_masked.gbff";   
    # launch snippy
    #/home/t-iris-005/SOFTWARE/snippy/bin/snippy --cpus 8 --outdir $OUTDIR --reference $REF --R1 $SAMPLE --R2 ${SAMPLE/_R1/_R2} ;
    /home/t-iris-005/SOFTWARE/snippy/bin/snippy --cpus 8 --outdir $OUTDIR --reference $REF --se $SAMPLE;
    #/home/t-iris-005/SOFTWARE/snippy/bin/snippy --cpus 8 --outdir $OUTDIR --reference $REF --R1 $SAMPLE --R2 ${SAMPLE/_R1/_R2} --mincov 3;
done


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
 ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
 ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

#-----------------------------------------------------------------------------------------
#################### STEP 02: Get Core SNP                        ########################
#-----------------------------------------------------------------------------------------

# step 0
# je propose d'utiliser le fichier xlsx avec une colonne panel pour choisir le panel que l'on souhaite!



# !! you need to be in the folder you want the output !!

# variables
SNIPPY_FOLDERS_LIST=$(cat list_44_strains.txt)  # will give all sample folders previously done
SNIPPY_FOLDERS_LIST=$(ls -d ../1-SNIPPY/*)
REF="/home/t-iris-005/0-RAW_DATA/References/Agy99_pMUM001_Hard_masked.gbff"
PREFIX="44_strains"

# launch snippy-core 3.2 (provide alignement but stats are bad!)
/home/t-iris-005/SOFTWARE/snippy-3.2/bin/snippy-core --prefix $PREFIX  $SNIPPY_FOLDERS_LIST



# launch snippy 4.3.8 (May 2019) don't work :( but provide nice stats like number of SNP to ref, missing snp...!
/home/t-iris-005/SOFTWARE/snippy/bin/snippy-core --prefix $PREFIX  $SNIPPY_FOLDERS_LIST --ref $REF



##### some supplementary stats on snippy-core results
# get SNP name
for tab in *tab; do cat $tab | cut -f 1,2 | sed 's/\t/_/g' | sed '1d' > ${tab}.SNP; done
# compare SNP name between two files (nombre de SNP commun)
comm -12 <(sort 179_strains_no_outgroup.tab.SNP) <(sort 202_strains_no_outgroup.tab.SNP) | wc -l


