################ Alban BESNARD ########### albanbesnard@hotmail.fr

#########################################################################################
######## Script shell pour traiter des données de séquençage         ####################
######## issues d'un méthode de capture innovante: Gexplore          ####################
######## GEXPLORE-MYCO-2019                                          ####################
#########################################################################################

#########################################################################################
###                                   SUMMARY                                         ###
#########################################################################################
###################
###################  00. Raw DATA
###################		0.1 Vérification de la qualité des données avec FASTQC
###################
###################  01. Filtering on quality and removing adaptators
###################
###################  02. MAPPING avec BWA on reference genome Agy99
###################
###################  03. PEAK_CALLING avec MACS
###################
##########################################################################################


#-----------------------------------------------------------------------------------------
#################### STEP 00: Raw data ########################
#-----------------------------------------------------------------------------------------

# On va dans le working directory
cd /home/t-iris-005/GEXPLORE-MYCO-2019

# On descend dans le dossier 00-RAW_DATA 
cd 00-RAW_DATA

###### 0.1 On vérifie la qualité des données avec un FASTQC
mkdir FASTQC
for i in *.gz; do
/home/t-iris-005/SOFTWARE/FastQC/fastqc $i -o FASTQC
done



## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
 ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
 ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

#-----------------------------------------------------------------------------------------
#################### STEP 01: cleaning data                     ########################
#-----------------------------------------------------------------------------------------

# on va nettoyer les données avec trimgalore
# quality filtre = 25.
# cela permet de supprimer les adaptateurs illumina et d'enlever
# les bases/reads de mauvaises qualités
cd ../01-CLEAN_DATA/

for R1 in ../00-RAW_DATA/*_R1_001.fastq.gz; do
R2=$(echo $R1 | sed 's/R1/R2/g')
/home/t-iris-005/SOFTWARE/TrimGalore-0.6.1/trim_galore --paired $R1 $R2 -q 25;
done



## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
 ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
 ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##


#-----------------------------------------------------------------------------------------
#################### STEP 02: mapping data on ref                     ####################
#-----------------------------------------------------------------------------------------
cd ../02-MAPPING/

# le mapping a proprement parlé
for i in ../01-CLEAN_DATA/*R1*; do
echo $i
indiv=$(basename $i | cut -f 2 -d _)
echo $indiv
bwa mem /home/t-iris-005/0-RAW_DATA/References/Agy99.fa $i ${i/R1_001_val_1/R2_001_val_2} > ${indiv}.sam
done

# condenser les .sam en .bam
SAM="S8.sam"
indiv=$(basename $SAM .sam)
BAM=${indiv}.bam
samtools view -Sb $SAM  > $BAM
samtools sort ${indiv}.bam -o ${indiv}.sorted.bam
samtools index ${indiv}.sorted.bam

rm ${indiv}.sam ${indiv}.bam


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
 ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
 ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##


#-----------------------------------------------------------------------------------------
#################### STEP 03: PEAK calling with Macs2                 ####################
#-----------------------------------------------------------------------------------------


mkdir ../03-PEAK_CALLING/
cd ../03-PEAK_CALLING/

MYCO_BAM=$(ls ../02-MAPPING/*.bam | head -4)
CONTROL_BAM=$(ls ../02-MAPPING/*.bam | tail -4)

# ICI on fait avec tous les bam dans le bon sens
macs2 callpeak -t $MYCO_BAM -c $CONTROL_BAM -f BAM -g ce -n all_bam -B -q 0.01

# ICI on fait avec tous les bam dans le bon sens
macs2 callpeak -t $MYCO_BAM -c $CONTROL_BAM -f BAM -g ce -n no_model -B -q 0.01 --nomodel --extsize 100

# ICI on fait avec tous les bam dans le bon sen
macs2 callpeak -t $CONTROL_BAM -c $MYCO_BAM -f BAM -g ce -n back_bam -B -q 0.01

# La on fait les bam myco seuls
macs2 callpeak -t $MYCO_BAM -f BAM -g ce -n myco_bam -B -q 0.01

# La on fait les bam control seuls
macs2 callpeak -t $CONTROL_BAM -f BAM -g ce -n control_bam -B -q 0.01

# chaque bam myco contre l'ensemble des témoins
for BAM in $MYCO_BAM; do 
echo $BAM;
indiv=$(basename $BAM | cut -f 1 -d _)
macs2 callpeak -t $BAM -c $CONTROL_BAM -f BAM -g ce -n $indiv -B -q 0.01
done












































