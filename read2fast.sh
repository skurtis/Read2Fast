#!/bin/bash

Ref=$1
FQ=$2
Sample=$3
BQ=$4
MQ=$5
MAC=$6
Mincov=$7

module load bwa/0.7.17 samtools/1.9

FILE=${Ref}1.amb
if [ ! -f $FILE ]; then
# If the reference FASTA hasn't been indexed and processed yet

awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $Ref | sed "/^$/d" > ${Ref}_intermediate.fasta
# The above command un-wraps the FASTA file so the whole sequence is on one line and removes spaces between lines.

wordcount="$(wc -l ${Ref}_intermediate.fasta | cut -d" " -f1)"
> ${Ref}1
i=1
while [ $i -le $wordcount ]
do
head -n $i ${Ref}_intermediate.fasta | tail -n1 > ${Ref}_current_line.txt
ifhead="$(grep -cE ">" ${Ref}_current_line.txt)"
if (( ifhead ))
then
cat ${Ref}_current_line.txt >> ${Ref}1
# If the current line is a header, forward it directly to the intermediate file.
else
cat ${Ref}_current_line.txt | sed -e 's/\?/N/g' | sed -e 's/-//g' >> ${Ref}1
# The above command replaces question marks with N's and removes gaps marked by dashes.
fi
i="$(( $i + 1 ))"
done

rm ${Ref}_current_line.txt ${Ref}_intermediate.fasta

bwa index -a is ${Ref}1
# Create a BWA index for the reference markers

fi

FILE2=${Ref}1.fai

if [ ! -f $FILE2 ]; then

samtools faidx ${Ref}1
# Make a .fai index file for the reference.

fi

fileno="$(ls $FQ | wc -l)"
fileno1="$(( $fileno - 1 ))"
if (( $fileno1 ))
then
bwa mem -t $SLURM_TASKS_PER_NODE -R '@RG\tID:'"$Sample"'\tSM:'"$Sample"'\tPL:ILLUMINA' ${Ref}1 ${FQ}/* | samtools sort -@ $SLURM_TASKS_PER_NODE -O bam -o ${Sample}_init.bam
# Reads are in two separate read files (reads 1 and 2)
else
bwa mem -t $SLURM_TASKS_PER_NODE -p -R '@RG\tID:'"$Sample"'\tSM:'"$Sample"'\tPL:ILLUMINA' ${Ref}1 ${FQ}/* | samtools sort -@ $SLURM_TASKS_PER_NODE -O bam -o ${Sample}_init.bam
# Reads are interleaved in one file
fi
# Perform BWA read-mapping of the sample reads to the reference markers.

samtools view -@ $SLURM_TASKS_PER_NODE -o ${Sample}_MQ.bam -h -b -q $MQ ${Sample}_init.bam 2>&1
# Select only reads that have a minimum MQ of 60

samtools rmdup ${Sample}_MQ.bam ${Sample}_MQrm.bam
# Remove PCR duplicates

rm ${Sample}_MQ.bam

module load gcc/5.2.0 freebayes/1.1.0-20170823

freebayes -b ${Sample}_MQrm.bam -f ${Ref}1 -v ${Sample}.vcf -n 2 --theta 0.001 --ploidy 2 --reference-quality '100,60' --haplotype-length 3 --min-repeat-size 5 --min-repeat-entropy 0 -m 1 -q $BQ -R 0 -Y 0 -e 1000 -F 0.0 -C $MAC -G 1 -Q 10 -U 1000 -z 1.0 --read-snp-limit 1000 --min-coverage $Mincov --min-alternate-qsum 0 --base-quality-cap 0 --prob-contamination 1e-08 -B '1000' -W '1,3' -D '0.9' --genotyping-max-banddepth 6

FILE3=${Ref}2

if [ ! -f $FILE3 ]; then
# If the reference FASTA hasn't been indexed and processed yet

wordcount="$(wc -l ${Ref}1 | cut -d" " -f1)"
i=1
> ${Ref}2
while [ $i -le $wordcount ]
do
head -n $i ${Ref}1 | tail -n1 > ${Ref}1_current_line.txt
ifhead="$(grep -cE ">" ${Ref}1_current_line.txt)"
if (( ifhead ))
then
cat ${Ref}1_current_line.txt >> ${Ref}2
# If the current line is a header, it is directly forwarded to the new file.
else
cat ${Ref}1_current_line.txt | sed -e 's/K\|W\|Y\|S\|R\|M\|B\|V\|D\|H/N/g' >> ${Ref}2
# The above command replaces any ambiguities with N's
fi
i="$(( $i + 1 ))"
done

rm ${Ref}1_current_line.txt

fi

module load gcc/5.2.0 bcftools/1.5 bamtools/2.1.1

# Remove unnecessary genotype information to prevent error
vcftools --vcf ${Sample}.vcf --recode --recode-INFO GT --out ${Sample}.vcf_GT 2>/dev/null

# Prevent reference alleles from being incorporated in the consensus FASTA
sed -ie 's/0\/1/1\/1/g' ${Sample}.vcf_GT.recode.vcf
sed -ie 's/0\/2/2\/2/g' ${Sample}.vcf_GT.recode.vcf
sed -ie 's/0\/[3-9]/0\/0/g' ${Sample}.vcf_GT.recode.vcf

# Generate a file of depth per each nucleotide.
samtools depth -a ${Sample}_MQrm.bam > ${Sample}_coverage.txt

# Create a tab-delimited file of each zero-coverage locus.
grep -E "[[:space:]]0" ${Sample}_coverage.txt > ${Sample}_zerocoverage.txt

cat ${Sample}_zerocoverage.txt | cut -f1,2 > ${Sample}_zerocov.txt

# BG-zip and index the VCF file
bgzip ${Sample}.vcf_GT.recode.vcf
bcftools index -c ${Sample}.vcf_GT.recode.vcf.gz

# Create a consensus FASTA file for each haplotype while masking zero-coverage loci.
bcftools consensus -H 1 -f ${Ref}2 -m ${Sample}_zerocov.txt ${Sample}.vcf_GT.recode.vcf.gz > ${Sample}_hap1.fasta
bcftools consensus -H 2 -f ${Ref}2 -m ${Sample}_zerocov.txt ${Sample}.vcf_GT.recode.vcf.gz > ${Sample}_hap2.fasta

rm ${Sample}.vcf_GT* ${Sample}_zerocoverage.txt ${Sample}_coverage.txt ${Sample}_zerocov.txt 

sed -i s/\>/\>${Sample}_hap1_/g ${Sample}_hap1.fasta
sed -i s/\>/\>${Sample}_hap2_/g ${Sample}_hap2.fasta
# Add sample name and haplotype to gene name

module load mafft

> ${Sample}_consensus.fasta
> ${Sample}_gaps.txt

fastalines="$(grep -c ">" ${Sample}_hap1.fasta)"
# Total number of sequences in FASTA file
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' ${Sample}_hap1.fasta | sed "/^$/d" > ${Sample}_hap1_unwrapped.fasta
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' ${Sample}_hap2.fasta | sed "/^$/d" > ${Sample}_hap2_unwrapped.fasta

j=1
while [ $j -le $fastalines ]
do

> ${Sample}_nucs.txt
head -n $(( $j * 2 )) ${Sample}_hap1_unwrapped.fasta | tail -n2 > ${Sample}_mafft.fasta
# Add haplotype 1 sequence to FASTA file to be aligned
head -n $(( $j * 2 )) ${Sample}_hap2_unwrapped.fasta | tail -n2 >> ${Sample}_mafft.fasta
# Add haplotype 2 sequence to FASTA file to be aligned

mafft --maxiterate 3 --globalpair ${Sample}_mafft.fasta > ${Sample}_mafft_global.fasta
# Use Mafft aligner to align the two haplotype sequences for that locus

awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' ${Sample}_mafft_global.fasta | sed "/^$/d" > ${Sample}_mafft_unwrapped.fasta
# Unwrap the Mafft alignment file
head -n1 ${Sample}_mafft_unwrapped.fasta | sed -e 's/_hap1//g' >> ${Sample}_consensus.fasta
head -n1 ${Sample}_mafft_unwrapped.fasta | sed -e 's/_hap1//g' >> ${Sample}_gaps.txt
# Output the header containing the sample and gene name into the growing consensus FASTA
grep -v ">" ${Sample}_mafft_unwrapped.fasta > ${Sample}_mafft_seqsonly.txt
# Grab only the sequences from the Mafft alignment file
grep -oE "\-{1,}" ${Sample}_mafft_seqsonly.txt | wc -l >> ${Sample}_gaps.txt
# Output the number of gaps present in alignment
chct="$(wc -c ${Sample}_mafft_seqsonly.txt | cut -d' ' -f1)"
chct2="$(( $chct / 2 - 1 ))"
# Total number of nucleotides

posn=1
while [ $posn -le $chct2 ]; do
cut -c $posn ${Sample}_mafft_seqsonly.txt > ${Sample}_posnucs.txt
# Select position in alignment
hd="$(head -n1 ${Sample}_posnucs.txt)"
# Grab haplotype 1 nucleotide
tl="$(tail -n1 ${Sample}_posnucs.txt)"
# Grab haplotype 2 nucleotide

if [[ $hd = $tl ]]; then
	echo $hd | tr '[:lower:]' '[:upper:]' >> ${Sample}_nucs.txt
# If the site is homozygous
else
# If the site is heterozygous or unresolved
        if [[ $hd = "t" ]]; then
                if [[ $tl = "c" ]]; then
                        echo "Y" >> ${Sample}_nucs.txt
                elif [[ $tl = "g" ]]; then
                        echo "K" >> ${Sample}_nucs.txt
                elif [[ $tl = "a" ]]; then
                        echo "W" >> ${Sample}_nucs.txt
                else
                        echo "N" >> ${Sample}_nucs.txt
                        # For all sites labeled as gaps (-) or unknown (N)
                fi
        elif [[ $hd = "a" ]]; then
                if [[ $tl = "t" ]]; then
                        echo "W" >> ${Sample}_nucs.txt
                elif [[ $tl = "c" ]]; then
                        echo "M" >> ${Sample}_nucs.txt
                elif [[ $tl = "g" ]]; then
                        echo "R" >> ${Sample}_nucs.txt
                else
                        echo "N" >> ${Sample}_nucs.txt
                fi
        elif [[ $hd = "c" ]]; then
                if [[ $tl = "t" ]]; then
                        echo "Y" >> ${Sample}_nucs.txt
                elif [[ $tl = "g" ]]; then
                        echo "S" >> ${Sample}_nucs.txt
                elif [[ $tl = "a" ]]; then
                        echo "M" >> ${Sample}_nucs.txt
                else
                        echo "N" >> ${Sample}_nucs.txt
                fi
        elif [[ $hd = "g" ]]; then
                if [[ $tl = "t" ]]; then
                        echo "K" >> ${Sample}_nucs.txt
                elif [[ $tl = "c" ]]; then
                        echo "S" >> ${Sample}_nucs.txt
                elif [[ $tl = "a" ]]; then
                        echo "R" >> ${Sample}_nucs.txt
                else
                        echo "N" >> ${Sample}_nucs.txt
                fi
        else
                echo "N" >> ${Sample}_nucs.txt
                # For all sites labeled as gaps (-) or unknown (N)
        fi
fi
posn="$(( $posn + 1 ))"
done

paste -s ${Sample}_nucs.txt | sed -e 's/\t//g' >> ${Sample}_consensus.fasta
# Align single nucleotides and add them to growing consensus FASTA

j="$(( $j + 1 ))"
done

awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' ${Sample}_consensus.fasta | sed "/^$/d" > ${Sample}_unwrappedconsensus.fasta
mv ${Sample}_hap1_unwrapped.fasta ${Sample}_hap1.fasta
mv ${Sample}_hap2_unwrapped.fasta ${Sample}_hap2.fasta
mv ${Sample}_unwrappedconsensus.fasta ${Sample}_consensus.fasta
rm ${Sample}_nucs.txt ${Sample}_mafft.fasta ${Sample}_mafft_seqsonly.txt ${Sample}_mafft_unwrapped.fasta ${Sample}_mafft_global.fasta ${Sample}_posnucs.txt 
