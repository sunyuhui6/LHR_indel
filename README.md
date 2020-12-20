# LHR_indel: LFR Haplotype-resolved Indel caller

LHR_indel is writen for individual genome phasing and assembly based on the stLFR and NGS technology. Latest version was v1.0 upload in 12/20/2020.

Usage:
perl Main.pl <fa.fai> <indexed sorted bam> <vcf> <window size> <max LFR> <min SupportBarcodes> <min Link> <shellDir> <outDir_tmp> <outDir_final> <fq1> <fq2>  <outDir_reads> <outDir_asm>  <refGenome>

Detailed Description:
<fa.fai> is the length file of the reference genome.
<indexed sorted bam> is an alignment file in BAM format which must be sorted and indexed in advance
<vcf> shoud store variant calling result in VCF format
<window size> is set for each phasing window, 10MB is suggested
<max LFR> is decided by the maximum length of DNA fragment. Normally it is not larger than 300KB according to the curren
stretage of stLFR technology.
<min SupportBarcodes> is the minimum count of supporting barcode for each variant
<min Link> is the minimum count of linked info
<shellDir> is an directory storing shell scripts which need to be executed, absolute path is suggested
<outDir_tmp> is a temporary directory storing split phasing result as well as barcode supporting info, absolute path is suggested
<outDir_final> is a directory storing the final phasing result in chr*.vcf files, absolute path is suggested.
<fq1> is the first fastq file 
<fq2> is the second fastq file 
<outDir_reads> is a directory storing clustered reads
<outDir_asm> is a directory storing final assembly contigs
<refGenome> is the reference genome 

Example:
perl Main.pl ref.fa.fai input.bam input.vcf 10000000 300000 1 1 shellDir outDir_tmp outDir_final fq1 fq2 outDir_reads outDir_asm refGenome
After running this command, you will find multiple shell scripts in shellDir, with the names of "run.chr1.sh" or "run.chrX.sh" or "run.scaffold1.sh"
which could be executed parallelly. In average, each chromosome costs 18~24 hours by consuming 4GB memory.The final phasing
results can be found in outDir_final dir in chr*.vcf files, while the medium files storing barcode supporting info can be
found in the outDir_tmp dir. The final assembly results can be fould in the outDir_asm dir.

Author:
Yuhui Sun

Email:
sunyuhui@genomics.cn
