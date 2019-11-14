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
ln -s /home/t-iris-005/0-RAW_DATA/B1815/*.fastq.gz  .
ln -s /home/t-iris-005/0-RAW_DATA/B1816/*.fastq.gz  .

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

for R1 in ../0-RAW_DATA/*B1816_R1_001.fastq.gz; do
R2=$(echo $R1 | sed 's/R1/R2/g')
/home/t-iris-005/SOFTWARE/TrimGalore-0.6.1/trim_galore --paired $R1 $R2 -q 25;
done

#### 1.2 checking quality of the data with FastQC
##########################################################
mkdir FASTQC
for i in *1.fq.gz; do
/home/t-iris-005/SOFTWARE/FastQC/fastqc $i -o FASTQC
donei

#-----------------------------------------------------------------------------------------
#################### STEP 2: mapping on reference ########################
#-----------------------------------------------------------------------------------------

#### 2.1 mapping using bwa
##########################################################

mkdir ../2-MAPPING
cd ../2-MAPPING/

# le mapping a proprement parlé
for i in ../../1-CLEAN_DATA/9*_R1_001_val_1*; do
echo $i
indiv=$(basename $i | cut -f 1,2 -d _)
echo $indiv
#bwa mem /home/t-iris-005/0-RAW_DATA/References/Agy99.fa $i ${i/R1_001_val_1/R2_001_val_2} > ${indiv}.sam
#bwa mem /home/t-iris-005/0-RAW_DATA/References/liflandii_128FXT.fa $i ${i/R1_001_val_1/R2_001_val_2} > ${indiv}.sam
#bwa mem /home/t-iris-005/0-RAW_DATA/References/ATCC33728.fasta $i ${i/R1_001_val_1/R2_001_val_2} > ${indiv}.sam
#bwa mem /home/t-iris-005/0-RAW_DATA/References/PlasmideCoreGenes.fasta $i ${i/R1_001_val_1/R2_001_val_2} > ${indiv}.sam
bwa mem /home/t-iris-005/0-RAW_DATA/References/IS2404-IS2606.fasta $i ${i/R1_001_val_1/R2_001_val_2} > ${indiv}.sam

# convert into .bam (less space used)
BAM=${indiv}.bam
samtools view -Sb ${indiv}.sam  > $BAM
samtools sort ${indiv}.bam -o ${indiv}.sorted.bam
samtools index ${indiv}.sorted.bam

rm ${indiv}.sam ${indiv}.bam

done

## merging all data from envrionnemental test (first and second run // WGA not WGA = 4 échantillons)
samtools merge merged_15-16.bam 15_B1816.sorted.bam 15.sorted.bam 16_B1816.sorted.bam 16.sorted.bam
samtools sort merged_15-16.bam -o merged_15-16.sorted.bam
samtools index merged_15-16.sorted.bam
rm merged_15-16.bam

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
#################### STEP 3: SNP detection ########################
#-----------------------------------------------------------------------------------------

cd ../3-SNP_DETECTION

# exporting varaibles and adding emboss to path
export PERL5LIB=/home/t-iris-005/.cpan/build/BioPerl-1.007002-A8KCCK
export PATH=$PATH:/home/t-iris-005/SOFTWARE/vt::/home/t-iris-005/Telechargements/EMBOSS-6.6.0/emboss::/home/t-iris-005/SOFTWARE/minimap2/

## for sample 7
SAMPLE="../1-CLEAN_DATA/7_S7_R1_001_val_1.fq.gz" 
OUTDIR="7"
REF="/home/t-iris-005/0-RAW_DATA/References/Agy99_pMUM001_Hard_masked.gbff"
# launch snippy
/home/t-iris-005/SOFTWARE/snippy/bin/snippy --cpus 8 --outdir $OUTDIR --reference $REF --R1 $SAMPLE --R2 ${SAMPLE/_R1_001_val_1/_R2_001_val_2} 


## for sample 8
SAMPLE="../1-CLEAN_DATA/8_S8_R1_001_val_1.fq.gz" 
OUTDIR="8"
REF="/home/t-iris-005/0-RAW_DATA/References/Agy99_pMUM001_Hard_masked.gbff"
# launch snippy
/home/t-iris-005/SOFTWARE/snippy/bin/snippy --cpus 8 --outdir $OUTDIR --reference $REF --R1 $SAMPLE --R2 ${SAMPLE/_R1_001_val_1/_R2_001_val_2} 

## for sample 13
SAMPLE="../1-CLEAN_DATA/13_S13_R1_001_val_1.fq.gz" 
OUTDIR="13"
REF="/home/t-iris-005/0-RAW_DATA/References/Agy99_pMUM001_Hard_masked.gbff"
# launch snippy
/home/t-iris-005/SOFTWARE/snippy/bin/snippy --cpus 8 --outdir $OUTDIR --reference $REF --R1 $SAMPLE --R2 ${SAMPLE/_R1_001_val_1/_R2_001_val_2} 



#-----------------------------------------------------------------------------------------
#################### STEP B: Kraken Krona ########################
#-----------------------------------------------------------------------------------------

# How many reads identified as anything?
# How many reads identified as mycobacteria?
# How many reads identified as Mycobacterium ulcerans?


#-----------------------------------------------------------------------------------------
#################### STEP C: IS2404 and IS2606 divesity in Agy99 ########################
#-----------------------------------------------------------------------------------------

cd C-IS2426_DIVERSITY
more Agy99.gff3 | grep "sequence:IS2404" | cut -f 1,4,5,7 | sed 's/Chromosome/CP000325.1/g' | sed 's/BX649209/BX649209.1/g' | sed 's/+/forward\t1\t+/' | sed 's/-/reverse\t1\t-/'> IS2404.bed
fastaFromBed -fi Agy99.fasta -bed IS2404.bed  -s -fo IS2404.fasta
more Agy99.gff3 | grep "sequence:IS2606" | cut -f 1,4,5,7 | sed 's/Chromosome/CP000325.1/g' | sed 's/BX649209/BX649209.1/g' | sed 's/+/forward\t1\t+/' | sed 's/-/reverse\t1\t-/'> IS2606.bed
fastaFromBed -fi Agy99.fasta -bed IS2606.bed  -s -fo IS2606.fasta
# alignement made using seaview "clustalo"
# choosing two Sequences to serve as model


## adding 4 choosen samples from diversity (2*MuA2 + 1*G1 +1*G5 + 2*G8 + 2*Cameroun)
echo "/home/t-iris-005/0-RAW_DATA/Genomes_Papous_Australie/282-5-11_R1.fastq.gz /home/t-iris-005/0-RAW_DATA/Genomes_Papous_Australie/320-7-12_R1.fastq.gz
/home/t-iris-005/0-RAW_DATA/Genomes_Cameroun/08-47_R1.fastq.gz /home/t-iris-005/0-RAW_DATA/Genomes_Cameroun/11-30_R1.fastq.gz
/home/t-iris-005/0-RAW_DATA/Genomes_Benin_samples/1357-12_R1.fastq.gz /home/t-iris-005/0-RAW_DATA/Genomes_Benin_samples/1448-13_R1.fastq.gz
/home/t-iris-005/0-RAW_DATA/Genomes_Benin_samples/1668-13_R1.fastq.gz /home/t-iris-005/0-RAW_DATA/Genomes_Benin_samples/1707-14_R1.fastq.gz" > list_diversity.txt

# mapping on reference
for i in $(cat list_diversity.txt); do
echo $i
indiv=$(basename $i | sed 's/_R1/;/' | cut -f 1 -d ";")
echo $indiv
bwa mem /home/t-iris-005/0-RAW_DATA/References/IS2404-IS2606.fasta $i ${i/R1_001_val_1/R2_001_val_2} > ${indiv}.sam

# convert into .bam (less space used)
samtools view -Sb ${indiv}.sam  > ${indiv}.bam
samtools sort ${indiv}.bam -o ${indiv}.sorted.bam
samtools index ${indiv}.sorted.bam

#mpilup (counting reads)
samtools mpileup ${indiv}.sorted.bam -A --reference /home/t-iris-005/0-RAW_DATA/References/IS2404-IS2606.fasta  > ${indiv}.mpileup
# getting variant only
varscan mpileup2cns ${indiv}.mpileup --output-vcf --min-coverage 8 --min-avg-qual 15 --min-var-freq 0.01 --min-reads2 2 --strand-filter 0 --p-value 1 --variants > ${indiv}.Varscan.vcf
# rename Sample
sed -i "s/Sample1/${indiv}/g" ${indiv}.Varscan.vcf
bgzip -c ${indiv}.Varscan.vcf > ${indiv}.Varscan.vcf.gz
tabix -p vcf ${indiv}.Varscan.vcf.gz




rm ${indiv}.sam ${indiv}.bam ${indiv}.mpileup ${indiv}.sorted.bam
done


# merging (not functionning!!) need to use the vcf file and the % of allele
bcftools merge *.vcf.gz > pool.Varscan.vcf
#rtg vcfsubset -i pool.Varscan.vcf -o pool_clean.Varscan.vcf --keep-format FREQ,PVAL -Z
rm pool_clean.Varscan.vcf
rtg vcfsubset -i pool.Varscan.vcf -o pool_clean.Varscan.vcf --keep-format FREQ -Z





