#!/bin/bash
#SBATCH --job-name=dada2do1       # Job name
#SBATCH --mail-user=<user@someadress.com> # Where to send mail
#SBATCH --mail-type=FAIL,END           	# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --partition=cpu-short		  	    # Partition to use for the analysis
#SBATCH --time=01:30:00               	# Time limit hrs:min:sec 
#SBATCH --ntasks=1				        # Number of tasks per node
#SBATCH --cpus-per-task=16				# Number of cpus requested per task
#SBATCH --export=none           		# Recommended for a clean working environment

# INPUT
script_dir=/home/Scripts/       # Location of all pipeline R scripts
fastq_dir=/data1/FastqFiles/    # Location of raw fastq files (forward files end with _1, reverse with _2)

Forward_Primer_Sequence=GTGARTCATCRARTYTTTG   # Known forward primer sequence
Reverse_Primer_Sequence=CCTSCSCTTANTDATATGC   # Known reverse primer sequence
Overwrite_Files=TRUE                          # Should previously generated files be overwritten (see manual for explanation)
Delete_Files=TRUE                             # Should N filtered reads be deleted afterwards to save space?

# DO NOT CHANGE
# Load modules and find path
module load R/3.6.2-fosscuda-2019b
module load cutadapt/1.18-GCCcore-8.2.0
cutadapt_path=`which cutadapt`

# Find script 
cd $script_dir

# Execute RScript 
time Rscript ./DADA2_Step_1.R --fastq_folder $fastq_dir --fw_primer $Forward_Primer_Sequence --rv_primer $Reverse_Primer_Sequence --overwrite $Overwrite_Files --delete $Delete_Files --cutadapt $cutadapt_path
