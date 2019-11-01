Read2Fast
=========
__This bash script maps sequencing reads to reference genomic loci, processes alignment files, produces the variant file, and generates consensus FASTAs.__

_Sarah Kurtis_

_Kimball/Braun Lab | Department of Biology | University of Florida_

_876 Newell Dr, Gainesville, FL 32611_

_[sarahkurtis@ufl.edu](mailto:sarahkurtis@ufl.edu)_

_September 19, 2019_
## Description
The following manual describes the steps and parameters of the program read2fast, which assembles trimmed sequencing reads for an input set of genomic loci and outputs BAM, VCF, and consensus haplotype and genotype FASTA files for latter genetic analyses. 
```
./read2fast.sh Ref FQ Sample BQ MQ MAC Mincov
```
Below is a description of the input parameters in order in which they must be entered.

__Parameter__ | __Description__
------------- | ---------------
__Ref__ | Reference FASTA file of genomic loci.
__FQ__ | Directory containing the trimmed reads FASTQ file(s), which may be either in the format of one file of paired reads 1 and 2 interleaved or two files of paired reads 1 and 2, respectively.
__Sample__ | Sample name to be used as the prefix for all output files and headers.
__BQ__ | Minimum base quality for variant calling. 
__MQ__ | Minimum base quality for variant calling. 
__MAC__ | Minimum alternate allele count per individual. 
__Mincov__ | Minimum coverage per nucleotide site. 
## Example
```
./read2fast.sh exons_reference.fasta trimmed_reads/ G_gallus_12 30 60 1 1
```
## Overview of steps 
1. The original reference FASTA file is run through a series of steps to ensure the following, which are requirements for the mapping and VCF-calling steps: 
    * No extra spaces between lines are allowed.
    * All gaps (denoted by dashes, or “-”) must be removed.
    * Unknown nucleotides denoted by “?” must be replaced with “N”. 
    * The FASTA must be un-wrapped such that all nucleotides are on one line per sequence.
2. The edited reference FASTA file has all IUPAC ambiguities converted to N’s, which is a requirement for the consensus-building step.
3. The first stage of the pipeline is the __bwa mem__ aligner set to be run on Illumina reads. The IS algorithm is used for indexing a set of genomic markers. The raw BAM file is outputted from the pipeline. 
4. The next step filters the initial BAM file using __samtools__ based on mapping quality of each alignment. The highly recommended MQ cutoff is 60. 
    * Due to the high rate of false positives that may be obtained when performing read mapping to reference markers especially of the same or closely-related species, it’s strongly suggested to download the assembly for a given marker and view a few sample alignments in Geneious or a similar software at varying mapping qualities. For highly-conserved loci, more than likely, only alignments with a MQ of 60 will appear as plausible alignments while others will be false positives. This initial filtering step significantly decreases the size of the file and makes subsequent analysis steps much quicker. 
5. The final cleanup step is to remove PCR duplicates. The cleaned BAM file is outputted from the pipeline.
6. The VCF file is generated from the final cleaned BAM file using __freebayes__ with several default parameters to note: 
    * The ploidy is set to 2.
    * The reference is assumed to not be of the same population as the sample.
    * All other parameters aside from those above and those specified as input to the pipeline are the default of __freebayes__.
7. __vcftools__ is then used to remove all unnecessary SNP information from the VCF file except for the genotype (GT) information to avoid latter errors. The VCF file is then BG-zipped and is indexed. 
8. As a requirement of the consensus builder, genotypes marked as “0/1” or “0/2” are converted to “1/1” or “2/2”, respectively, to prevent the reference allele from being outputted in the sample consensus FASTA. 
9. The pipeline generates a file of sequencing depth per nucleotide. Zero-coverage loci are extracted into a separate file. 
10. __bcftools consensus__ overlays the sample variants from the VCF file with the reference FASTA to generate a consensus FASTA file for each of the two haplotypes for the given sample, which is outputted from the pipeline. The masking option uses the zero-coverage loci file to convert any loci with zero coverage to N’s. The sequence header names in the consensus FASTAs are edited to match the sample prefix. 
11. The two haplotype FASTA files are then aligned using the MAFFT aligner and heterozygous loci are called. Insertions/deletions between haplotypes are called as N’s. This genotype FASTA file is outputted from the pipeline.  
