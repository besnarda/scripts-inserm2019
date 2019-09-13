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
########  0. Raw DATA
########     0.1 symbolic links with the data
########     0.2 checking data quality with fastqc 
########
########  1. cleaning adaptators and quality
########     1.1 cleaing with trimgalore
########     1.2 checking data quality with fastqc 
########
########  2. Mapping on ulcerans reference (Agy99)
########	 2.1 actual mapping using bwa
########     2.2 some basic statistic and visualisation
########		2.2.1 percentage of mapped reads
########        2.2.2 cover all over the genome
########
########  3. SNP detection (if possible)
########
########  B. Analysis of kraken krona done on galaxy (not enough memory on my computer)
##########################################################################################


#-----------------------------------------------------------------------------------------
#################### STEP 0: Raw data ########################
#-----------------------------------------------------------------------------------------

#### 0.1 importing data used for analysis
##########################################################

cd /home/t-iris-005/ENVIRO-MU-GENOME-2019
mkdir 0-RAW_DATA
cd 0-RAW_DATA
ln -s /home/t-iris-005/0-RAW_DATA/Genomes_Papous_Australie/282-5-11_R*.fastq.gz  .

#### 0.2 checking quality of the data with FastQC
##########################################################

mkdir FASTQC
for i in *1.fastq.gz; do
/home/t-iris-005/SOFTWARE/FastQC/fastqc $i -o FASTQC
done

#-----------------------------------------------------------------------------------------
#################### STEP 1: clean data ########################
#-----------------------------------------------------------------------------------------

#### 1.1 Cleaning data with Trimgalore
##########################################################

# on va nettoyer les données avec trimgalore
# quality filter = 25.

mkdir ../1-CLEAN_DATA
cd ../1-CLEAN_DATA/

for R1 in ../0-RAW_DATA/*_R1.fastq.gz; do
R2=$(echo $R1 | sed 's/R1/R2/g')
/home/t-iris-005/SOFTWARE/TrimGalore-0.6.1/trim_galore --paired $R1 $R2 -q 25;
done

#### 1.2 checking quality of the data with FastQC
##########################################################
mkdir FASTQC
for i in *1.fq.gz; do
/home/t-iris-005/SOFTWARE/FastQC/fastqc $i -o FASTQC
done

#-----------------------------------------------------------------------------------------
#################### STEP 2: mapping on reference ########################
#-----------------------------------------------------------------------------------------

#### 2.1 mapping using bwa
##########################################################

mkdir ../2-MAPPING
cd ../2-MAPPING/

# le mapping a proprement parlé
for i in ../1-CLEAN_DATA/*R1_*; do
echo $i
indiv=$(basename $i | cut -f 1 -d _)
echo $indiv
bwa mem /home/t-iris-005/0-RAW_DATA/References/Agy99.fa $i ${i/R1_val_1/R2_val_2} > ${indiv}.sam

# convert into .bam (less space used)
BAM=${indiv}.bam
samtools view -Sb ${indiv}.sam  > $BAM
samtools sort ${indiv}.bam -o ${indiv}.sorted.bam
samtools index ${indiv}.sorted.bam

rm ${indiv}.sam ${indiv}.bam
done

#### 2.2 some basic statistics
##########################################################


for i in ../1-CLEAN_DATA/*R1_*; do
	indiv=$(basename $i | cut -f 1 -d _)
	echo $indiv
	# nombre de reads au début (avant filtrage)
	a=$(more ../1-CLEAN_DATA/${indiv}_R1.fastq.gz_trimming_report.txt | grep "Total reads processed" | cut -f 2 -d ":" | sed 's/,//g' | sed 's/ //g')
	b=$(more ../1-CLEAN_DATA/${indiv}_R2.fastq.gz_trimming_report.txt | grep "Total reads processed" | cut -f 2 -d ":" | sed 's/,//g' | sed 's/ //g')
	nread_step0=$(($a+$b))
	echo $nread_step0

	# nombre de reads après filtrage
	a=$(more ../1-CLEAN_DATA/${indiv}_R1.fastq.gz_trimming_report.txt | grep "Reads written (passing filters)" | cut -f 2 -d ":" | cut -f 1 -d "(" | sed 's/,//g' | sed 's/ //g')
	b=$(more ../1-CLEAN_DATA/${indiv}_R2.fastq.gz_trimming_report.txt | grep "Reads written (passing filters)" | cut -f 2 -d ":" | cut -f 1 -d "(" | sed 's/,//g' | sed 's/ //g')
	nread_step1=$(($a+$b))
	echo $nread_step1
	# nombre de reads mappant sur ulcerans
	nread_step2=$(samtools flagstat ../2-MAPPING/${indiv}.sorted.bam | head -n 5 | tail -n 1 | cut -f 1 -d ' ')
	echo $nread_step2
done


#-----------------------------------------------------------------------------------------
#################### STEP B: Kraken Krona ########################
#-----------------------------------------------------------------------------------------

# How many reads identified as anything?
# How many reads identified as mycobacteria?
# How many reads identified as Mycobacterium ulcerans?