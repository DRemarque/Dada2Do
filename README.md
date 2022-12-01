# Dada2Do
A user friendly dada2 based pipeline for remote and local machines. Here, various methods and workflows have been collected into 
one automated pipeline that can do the following:
* Take minimised user input at the start and then run independantly
* Recognise previously made files to continue after unforseen stops
* Visualise data quality to help determine filtering 
* Use any quality filtering parameters 
* Automatically log used filtering criteria and how many reads they removed
* Deal with varying amplicon lengths of i.e. ITS regions
* Output ASVs or any percentage similarity OTUs
* Run pooled/unpooled samples (or both in parallel)
* Assign taxonomy with different databases in parallel

Without requiring the end user to know:
* More than basic R coding
* What base filtering should always be performed
* Different sequence error estimation approaches
* How to visualise data quality
* How to deal with varying amplicon sizes
* How to properly cluster ASVs to OTUs
* How to clean and manage abundance and taxonomy tables

So what does the pipeline require?
* Known primer sequences
* Access to cutadapt (https://cutadapt.readthedocs.io/en/stable/)
* Fastq files with _1 as forward and _2 as reverse (i.e. somename_1.fq.gz or somename_1.fastq, but not somename_fwread.fq.gz))
* One folder with all fastq files to be analysed together and access to write new files into that folder.

# How to use the pipeline?
Please view the manual in the documentation folder

# R package dependencies
BiocManager
dada2
optparse
dplyr
DECIPHER
Biostrings
msa
ggplot2
ggpubr
RcppParallel
openxlsx
ShortRead
seqTools

# Output files and steps
![image](https://user-images.githubusercontent.com/67581284/205045420-971d5a80-85bd-4264-8dae-1f406dafbc9d.png)


