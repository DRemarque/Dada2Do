################################################################################
#                           DADA2 Pipeline STEP 2                              #
################################################################################
message(
"
This is the second step in the DADA2 pipeline to go from single sample fastq files
with no primers/adapters to ASV and OTU tables.
Please ensure the R working directory is set to the R script locations!
")

#### Collect all the functions ####
source('DADA2_Step_2_Source_Code.R')

#### Collect the given arguments ####
arguments <- Argument_Parsing()
## Unsure what an argument is/does? Unhash and run below for help:
# Argument_Parsing(show_help=T)

#### When running locally, you can input your directories below and unhash ####

# arguments$fastq_folder <- '/home/User/Myfolder/'
# arguments$database     <- '/home/User/Mydatabase.fa.gz'

## All other arguments can be left unchanged if the defaults below are okay
## Otherwise change them as necessary and unhash

# arguments$novaseq      <- FALSE
# arguments$pool_samples <- FALSE
# arguments$overwrite    <- TRUE

# arguments$maxEE        <- 2
# arguments$truncLen_fw  <- 0
# arguments$truncLen_rv  <- 0
# arguments$minLen       <- 20
# arguments$truncQ       <- 2

# arguments$otu_similarity <- 1
# arguments$otu_abundance  <- TRUE

#### Run the code ####
main(opt=arguments)



