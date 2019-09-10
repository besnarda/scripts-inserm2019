################ Alban BESNARD ########### albanbesnard@hotmail.fr

#########################################################################################
######## Shell (=Linux language) script to redo all bioinformatic analysis ##############
######## from the environemental Whole genome sequencing of ulcerans host  ##############
######## Sequencing was done after a capture                               ##############
######## FOLDER NAME: ENVIRO_MU_GENOME_2019								   ##############
#########################################################################################

#########################################################################################
###                                   SUMMARY                                         ###
#########################################################################################
########
########  00. Raw DATA
########     0.1 checking data quality with fastqc
########     0.2 cleaning adaptators and quality
########
########  01. Mapping on ulcerans reference (Agy99)
########	 1.1 actual mapping using bowtie2
########     1.2 some basic statistic and visualisation
########		1.2.1 pourcentage de reads mappés sur ulcerans
########        1.2.2 graph couverture selon position
########     1.3 keeping reads which mapped for next steps
########
########  02. SNP detection (if possible)
##########################################################################################


#-----------------------------------------------------------------------------------------
#################### STEP 00: Raw data ########################
#-----------------------------------------------------------------------------------------

#### 0.1 checking quality of the date with FastQC
##########################################################

mkdir FASTQC
for i in *1.fastq.gz; do
/home/t-iris-005/SOFTWARE/FastQC/fastqc $i -o FASTQC
done

#### 0.2 Cleaning data with Trimgalore
##########################################################

# on va nettoyer les données avec trimgalore
# quality filter = 25.
# This is 
cd ../01-CLEAN_DATA/

for R1 in ../00-RAW_DATA/*_R1_001.fastq.gz; do
R2=$(echo $R1 | sed 's/R1/R2/g')
/home/t-iris-005/SOFTWARE/TrimGalore-0.6.1/trim_galore --paired $R1 $R2 -q 25;
done