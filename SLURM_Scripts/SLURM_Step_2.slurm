#!/bin/bash
#SBATCH --job-name=dada2do2               # Job name
#SBATCH --mail-user=<user@someadress.com>    # Where to send mail
#SBATCH --mail-type=FAIL,END           			# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --partition=cpu-medium                  # Partition to use for the analysis
#SBATCH --time=10:00:00                         # Time limit hrs:min:sec (you may not want this)
#SBATCH --ntasks=1                              # Number of tasks per node
#SBATCH --cpus-per-task=12                      # Number of cpus requested per task
#SBATCH --export=none                           # Recommended for a clean working environment

# INPUT
script_dir=/home/Scripts/               # Location of all pipeline R scripts
fastq_dir=/data1/Fastq_Files/           # Location of raw fastq files
database=/home/database.fa.gz           # Location of reference database (can be .gz)

Is_Novaseq=TRUE              # The sequencing method was Novaseq (TRUE) or something else like MiSeq (FALSE)

maxEE=2                      # Maximum estimated error in a read before it is discarded
truncLen_forward_reads=0     # Maximum bp forward read lenght before truncation. Shorter reads are discarded
truncLen_reverse_reads=0     # Maximum bp reverse read lenght before truncation. Shorter reads are discarded
minLen_read=20               # When truncLen=0 (turned off). How long should the minimal read bp be
truncQ=2                     # When truncLen=0 (turned off). Truncate reads when they drop below this quality score

Pool_Samples=FALSE           # When identifying ASVs, should samples be pooled? TRUE = high rare sequence sensitivity, FALSE = high between sample differentiation, PSEUDO = middle ground

otu=1                        # Percentage similarity for OTU's. 1 = no OTU, 0.97 = 97% OTU
otu_abundance=TRUE           # Whether the OTU representative sequence is most abundant (TRUE) or consensus (FALSE)

Overwrite_Files=FALSE        # Should previously generated files be overwritten (see manual for explanation)


# DO NOT CHANGE
# Load modules
module load R/3.6.2-fosscuda-2019b

# Find script 
cd $script_dir

# Execute RScript 
time Rscript ./DADA2_Step_2.R --fastq_folder $fastq_dir --database $database --novaseq $Is_Novaseq --maxEE $maxEE --truncLen_fw $truncLen_forward_reads --truncLen_rv $truncLen_reverse_reads --minLen $minLen_read --truncQ $truncQ --otu_similarity $otu --otu_abundance $otu_abundance --pool_samples $Pool_Samples --overwrite $Overwrite_Files
