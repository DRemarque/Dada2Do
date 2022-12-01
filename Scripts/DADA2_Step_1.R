################################################################################
#                           DADA2 Pipeline STEP 1                              #
################################################################################
message(
"
This is the first step in the DADA2 pipeline to go remove primer/adapter sequences and
to explore FASTQ file quality. 
Please ensure the R working directory is set to the R script locations!
")

#### Collect all the functions ####
source('DADA2_Step_1_Source_Code.R')

#### Collect the given arguments ####
arguments <- Argument_Parsing()
## Unsure what an argument is/does? Unhash and run below for help:
# Argument_Parsing(show_help=T)

#### When running locally, you can input your directories below and unhash ####

# arguments$fastq_folder <- 'C:/Users/drema/Documents/1.Uni/10EC_Proj/Test_Data'
# arguments$cutadapt     <- 'C:/Path/to/cutadapt.exe'
# arguments$fw_primer    <- 'GTGARTCATCRARTYTTTG'
# arguments$rv_primer    <- 'CCTSCSCTTANTDATATGC'

## All other arguments can be left unchanged if the defaults below are okay
## Otherwise change them as necessary and unhash

# arguments$overwrite    <- TRUE
# arguments$delete       <- TRUE

#### Run the code ####
main(opt=arguments)


