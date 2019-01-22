# Trimmitomics

#Data acquisition from NCBI Nucleotide Archive in SRA file format 
The datasets were converted to FASTQ using fastq-dump (with –readids and –split-files parameters) from SRAtoolkit.2.8.2.

#check the quality of your reads with FastQC; if you have more files use MultiQC
fastqc file_R1.fq.gz

#clean the reads
java -jar ~/trimmomatic-0.33.jar PE -threads 8 -phred33 ./file_R1.fq.gz ./file_R2.fq.gz file_R1_paired.fq.gz file_R1_unpaired.fq.gz file_R2_paired.fq.gz file_R2_unpaired.fq.gz ILLUMINACLIP:Adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:30&

#make a bowtie_index folder
bowtie2-build /path/to/file/REFgenome.fa bowtie_index/nameit

#alignment either --local or --end-to-end
bowtie2 --fr --local -x bowtie_index/nameit -1 file_R1_paired.fq.gz -2 file_R2_paired.fq.gz -S filename_BOWTIE.sam&

#Convert SAM to BAM using samtools 
samtools view -b -S -o filename.bam filename_BOWTIE.sam&

#extracting only the mapped reads 
samtools view -b -F 4 filename.bam > filename_mapped.bam&

#genome guided Trinity following instructions on the page https://github.com/trinityrnaseq/trinityrnaseq/wiki/Genome-Guided-Trinity-Transcriptome-Assembly
samtools sort filename_mapped.bam -o filename_mapped_sorted_Trinity.bam&

~/Trinity --genome_guided_bam /path/to/file/filename_Trinity.bam --max_memory 50G --genome_guided_max_intron 1000 --CPU 6 --output /path/to/file/trinity_out_dir --full_cleanup& > trinity.log&

#output is txt file with your contigs in *_GG.txt
