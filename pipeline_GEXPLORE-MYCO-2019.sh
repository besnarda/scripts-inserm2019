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

CONTROL_BAM=$(ls ../02-MAPPING/*.bam | head -4)
MYCO_BAM=$(ls ../02-MAPPING/*.bam | tail -4)

# ICI on fait avec tous les bam dans le bon sens
macs2 callpeak -t $MYCO_BAM -c $CONTROL_BAM -f BAM -g ce -n all_bam -B -q 0.01

# ICI on fait avec tous les bam dans le bon sens
macs2 callpeak -t $MYCO_BAM -c $CONTROL_BAM -f BAM -g ce -n size_200_keep_dup -B -q 0.01 --nomodel --extsize 200 --keep-dup all

# ICI on fait avec tous les bam dans le mauvais sens
macs2 callpeak -t $CONTROL_BAM -c $MYCO_BAM -f BAM -g ce -n back_bam -B -q 0.01 --nomodel --extsize 200 --keep-dup all

# La on fait les bam myco seuls
macs2 callpeak -t $MYCO_BAM -f BAM -g ce -n myco_bam -B -q 0.01 --nomodel --extsize 200 --keep-dup all

# La on fait les bam control seuls
macs2 callpeak -t $CONTROL_BAM -f BAM -g ce -n control_bam -B -q 0.01 --nomodel --extsize 200 --keep-dup all

# chaque bam myco contre l'ensemble des témoins
for BAM in $MYCO_BAM; do 
echo $BAM;
indiv=$(basename $BAM | cut -f 1 -d _)
macs2 callpeak -t $BAM -c $CONTROL_BAM -f BAM -g ce -n $indiv -B -q 0.01
done


### mycolactone VS temoin

input='size_200_keep_dup_peaks.xls'
>associated_gene.txt
while read LINE;do 
    echo $LINE
    chrom=$(echo $LINE | cut -f 1 -d " " )
    if [ $chrom = "Chromosome" ]; then chrom="CP000325"; fi
    mid=$(echo $LINE | cut -f 5 -d " " )
    tmp=$(cat ~/0-RAW_DATA/References/Agy99.gff3 | grep $chrom | grep -P '\tgene\t|\tpseudogene\t' | awk -v a="$mid" ' {if($4 < a && $5 > a) {print $4,$5,$7,$9}}' | sed 's/ /\t/g' | head -n1)
    echo $tmp >> associated_gene.txt
done < $input

paste size_200_keep_dup_peaks.xls associated_gene.txt |  sed 's/ /\t/g' > size_200_keep_dup_annotated.txt

### mycolactone

input='back_bam_peaks.xls'
>associated_gene.txt
while read LINE;do 
    chrom=$(echo $LINE | cut -f 1 -d " " )
    if [ $chrom = "Chromosome" ]; then chrom="CP000325"; fi
    mid=$(echo $LINE | cut -f 5 -d " " )
    tmp=$(cat ~/0-RAW_DATA/References/Agy99.gff3 | grep $chrom | grep -P '\tgene\t|\tpseudogene\t' | awk -v a="$mid" ' {if($4 < a && $5 > a) {print $4,$5,$7,$9}}' | sed 's/ /\t/g' | head -n1)
    echo $tmp >> associated_gene.txt
done < $input

paste back_bam_peaks.xls associated_gene.txt |  sed 's/ /\t/g' > back_bam_peak_annotated.txt


### témoin

input='control_bam_peaks.xls'
>associated_gene.txt
while read LINE;do 
    echo $LINE
    chrom=$(echo $LINE | cut -f 1 -d " " )
    if [ $chrom = "Chromosome" ]; then chrom="CP000325"; fi
    mid=$(echo $LINE | cut -f 5 -d " " )
    tmp=$(cat ~/0-RAW_DATA/References/Agy99.gff3 | grep $chrom | grep -P '\tgene\t|\tpseudogene\t' | awk -v a="$mid" ' {if($4 < a && $5 > a) {print $4,$5,$7,$9}}' | sed 's/ /\t/g' | head -n1)
    echo $tmp >> associated_gene.txt
done < $input

paste control_bam_peaks.xls associated_gene.txt |  sed 's/ /\t/g' > control_peak_annotated.txt

more myco_peak_annotated.txt | grep gene | cut -f 4 -d "-" | cut -f 1 -d ";"

####### passage du .xls au peptide d'intérêt.
name="back_bam"
# passage .xls au format .bed
cat ${name}_peaks.xls | grep -v "#" | cut -f 1,2,3 | tail -n +3 > ${name}_annotated.bed

# extraction des nucléotides (ne pas oublier de modifier le .xls en .txt format bed 3 colonnes tabulation & no header).
fastaFromBed -fi ~/0-RAW_DATA/References/Agy99.fa -bed ${name}_annotated.bed  -fo ${name}.fasta

# a small modification
cat ${name}.fasta | sed 's/:/_/g' > tmp.txt
mv tmp.txt ${name}.fasta 

# traduction
transeq ${name}.fasta -outseq ${name}_prot.fasta -frame 6

# get only prot with no stop codon
cat ${name}_prot.fasta | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | awk -F '\t'  '!($2 ~ /\*/)' | tr "\t" "\n" > ${name}_prot_filtered.fasta

# some stat on these peptides
pepstats ${name}_prot_filtered.fasta -outfile ${name}_prot_filtered.pepstats

# combien de peak ont une protéine associé sans codon stop?
cat ${name}_prot_filtered.fasta | grep ">" | cut -f 1,2 -d _ | uniq -c | wc -l

# combien de peptide potentiel?
cat ${name}_prot_filtered.fasta | grep -c ">" 



#-----------------------------------------------------------------------------------------
#################### STEP 05: Analyse des peaks sur mlsA1 et mlsB     ####################
#-----------------------------------------------------------------------------------------

# peak.bed est un fichier avec le nom du gène et début et fin. 
fastaFromBed -fi mls_prot.fasta -bed peak.bed  -fo mls_KS_AT.fasta

# on observe clairement une différence entre les modules qui matchent avec la myco et les autres!

#>mlsB:513-572
#HPHRATITTSIEHHSENNHDTTDALAALHALANNGTHPLLSRGLLTPQGPGKTVFVFPG
#>mlsB:2319-2365
#RAVVVGADRHQLQRGLAELASGNLGADVVVGRARAAGETVMVFPGQ






















