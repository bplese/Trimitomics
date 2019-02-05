# Trimitomics: An efficient pipeline for mitochondrial assembly from transcriptomic reads in non-model species
Author(s): [Bruna Plese], Maria Eleonora Rossi, Nathan James Kenny, Sergi Taboada, Vasiliki Koutsouveli, and Ana Riesgo (mailto:bplese@irb.hr)  

## Data

The datasets used are deposited at NCBI Sequence Read Archive (SRA) https://www.ncbi.nlm.nih.gov/sra. Conversion from SRA file format to FASTQ is done using fastq-dump (with –readids and –split-files parameters) from SRAtoolkit.2.8.2. All the reads used are paired end but single end reads can be used as well. Check the quality of reads using FastQC of for more files you can use MultiQC.

## Requisites

* NOVOPlasty version 2.7.1 
* Bowtie2 
* Trinity version 2013_08_14
* Velvet 1.2.10

## Mitochondrial genome recovery from RNA-seq data
* The proposed pipline uses 3 assemblies:

1. NOVOPlasty (https://github.com/ndierckx/NOVOPlasty)
2. Bowtie2/Trinity (https://github.com/trinityrnaseq/trinityrnaseq/wiki/Genome-Guided-Trinity-Transcriptome-Assembly)
3. Velvet (https://github.com/dzerbino/velvet)

## 1. NOVOPlasty 
* use suitable seed or reference genome for assembly 
* try different k-mer sizes (25, 39, 45, 51) in order to get better results
* output file contains either complete circular genome or contigs needed to be assembled
* in case you get contigs further assembly is required

## 2. Bowtie2/Trinity

Using Bowtie2 for alignement of raw reads with reference mitochondrial genome. Extracting reads that mapped to the reference mitochondrial genome. 

* clean the reads:
java -jar /path/trimmomatic-0.33.jar PE -threads 8 -phred33 /path/file_R1.fq.gz /path/file_R2.fq.gz /path/file_R1_paired.fq.gz /path/file_R1_unpaired.fq.gz /path/file_R2_paired.fq.gz /path/file_R2_unpaired.fq.gz ILLUMINACLIP:Adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:30&

* make a bowtie_index folder:
bowtie2-build /path/REFgenome.fa bowtie_index/name

* alignment --local with paired reads:
bowtie2 --fr --local -x bowtie_index/name -1 /path/file_R1_paired.fq.gz -2 /path/file_R2_paired.fq.gz -S /path/filename_BOWTIE.sam&

* Convert SAM to BAM using samtools: 
samtools view -b -S -o /path/filename.bam /path/filename_BOWTIE.sam&

* Extract only mapped reads: 
samtools view -b -F 4 /path/filename.bam > /path/filename_mapped.bam& #file size of output .bam file is small enough for easier further manipulation

* genome guided Trinity assembly:

samtools sort /path/filename_mapped.bam -o /path/filename_mapped_sorted.bam&

/path/Trinity --genome_guided_bam /path/filename_mapped_sorted.bam --max_memory 50G --genome_guided_max_intron 1000 --CPU 6 --output /path/trinity_out_dir --full_cleanup& > trinity.log&

*adjust intron size depending on the species with genome_guided_max_intron 1000
*mitochondrial contigs needed to be further assembled are in the folder trinity_out_dir; file name *_GG.txt


## 3. Velvet assembly

* run velvet for the beginning with a k-mer size 71: 
velveth velvet71 71 -shortPaired -separate -fastq.gz /path/file_R1_paired.fq.gz /path/file_R2_paired.fq.gz -short -fastq.gz /path/file_all_unpaired.fq.gz&

* this you use if you want to use only paired reads:
velveth velvet71 71 -shortPaired -separate -fastq.gz /path/file_R1_paired.fq.gz /path/file_R2_paired.fq.gz&

* if you would like to include unpaired as well: 
velveth velvet71 71 -shortPaired -separate -fastq.gz /path/file_R1_paired.fq.gz /path/file_R2_paired.fq.gz /path/file_all_unpaired.fq.gz&

* run velvetg, with a low coverage cutoff, just to exclude mistakes:
*if you work in the same folder there is no need for setting the path. 
velvetg velvet71/ -min_contig_lgth 100 -cov_cutoff 3 > jobname.log&

* In order to set a cutoff you blast the contigs that you got from velvetg above (in the velvet71 folder, file name contigs), using cox1 seq of the most closest species you find:

makeblastdb -in contigs.fa -title name_coxDB -out name_coxDB -dbtype nucl& 

blastn -query knownseq_cox1.fa -db name_coxDB -out nameofthefile.txt&

* or you can make blast with protein seq. we used this:

makeblastdb -in file_name_contigs.fa -title name_coxDB -out name_coxDB -dbtype nucl&

tblastn -query knownproteinseq_cox1.fa -db name_coxDB -db_gencode 4 -out nameofthefile.txt&

*-db_gencode 4 this is translation code for sponge mitochondria; adjust it according to the species you would like to assemble
*makeblastdb=make a database 

velvetg velvet71/ -min_contig_lgth 100 -cov_cutoff 80 -max_coverage 170&

*after that output file is contigs.fa in the folder velvet71
*In the final contigs.fa if you are lucky you will have mitodata in few contigs. If not, repeat step with blast, using known protein seq. of coding genes from the closest species possible. 


*if you still have gaps in your mitogenomes you can adjust the following:
1. -k-mer size 31,51,71
2. -cov_cutoff and -max_coverage 


## Final assembly
If the mt genome is not recovered by any of the three methods, the results are combined as a meta-assembly in order to obtain the best results. This can be done using any alignment software. 

## Checkpoint
When the complete or almost complete mt genome is obtained, the assembly data is checked for homology with BlastN against the NCBI nr database. Further annotation is then performed on the web server MITOS2 (Bernt et al., 2013) using the appropriate translation table.
