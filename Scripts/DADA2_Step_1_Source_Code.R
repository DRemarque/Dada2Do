message("
This script is used for multisample, demultiplexed sequence data with known 
primer sequences. \n
It checks and removes the primer sequences and helps determining data quality and filter thresholds \n
Please make sure each forward read file ends with _1, and each reverse file with _2 in the name.")

#### Collect all functions used in main ####
Argument_Parsing<-function(show_help=F){
  # Index Input
  option_list = list(
    optparse::make_option(c("--fastq_folder"), type="character", default=NULL,
                          help="The directory with all fastq files to be analysed",metavar = "character"),
    optparse::make_option(c("--fw_primer"), type="character", default=NULL,
                          help="Known forward primer used",metavar = "character"),
    optparse::make_option(c("--rv_primer"), type="character", default=NULL,
                          help="Known reverse primer used",metavar = "character"),
    optparse::make_option(c("--overwrite"), type="logical", default=T,
                          help="Whether previously made files should be overwritten",metavar = "logical"),
    optparse::make_option(c("--delete"), type="logical", default=T,
                          help="Wheter N filtered FASTQ files should be deleted to save disk space",metavar = "logical"),
    optparse::make_option(c("--cutadapt"), type="character", default=NULL,
                          help="Server path to cutadapt",metavar = "character")
    )
  opt_parser <- optparse::OptionParser(option_list=option_list)
  if(show_help){
    optparse::print_help(opt_parser)
  } else{
    opt <- optparse::parse_args(opt_parser)
    return(opt)  
  }
}

Install_Check<-function(){
  installed <- installed.packages()
  cran <- 'https://mirror.lyrahosting.com/CRAN/'
  if (!'BiocManager' %in% installed) {
    message('Bioconductor not found. Installing now\n')
    suppressMessages(install.packages("BiocManager",repos=cran))
  }
  if (!'dada2' %in% installed) {
    message('dada2 not found. Installing now\n')
    suppressMessages(BiocManager::install("dada2"))
  }
  if (!'ShortRead' %in% installed) {
    message('ShortRead not found. Installing now\n')
    suppressMessages(BiocManager::install("ShortRead"))
  }
  if (!'seqTools' %in% installed) {
    message('seqTools not found. Installing now\n')
    suppressMessages(BiocManager::install("seqTools"))
  }
  if (!'optparse' %in% installed) {
    message('optparse not found. Installing now\n')
    suppressMessages(install.packages('optparse',repos=cran))
  }
  
  message('All dependencies available\n')
}

allOrients <- function(primer) {
  # Create all orientations given primers
  dna <- Biostrings::DNAString(primer)  
  orients <- c(Forward = dna, 
               Complement = Biostrings::complement(dna), 
               Reverse = Biostrings::reverse(dna), 
               RevComp = Biostrings::reverseComplement(dna))
  return(sapply(orients, Biostrings::toString))  # Convert back to character vector
}

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- Biostrings::vcountPattern(primer, ShortRead::sread(ShortRead::readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

plot_read_length <- function(read_len,output_file=NA){
  
  m <- order(read_len$num_reads,decreasing=T)[1:5]
  l <- which(read_len$num_reads>10)[1]
  
  p1 <-ggplot2::ggplot(data=read_len,ggplot2::aes(x=.data$read_length,y=.data$num_reads)) + 
    ggplot2::geom_bar(stat="identity") +
    ggplot2::geom_text(data=read_len[m,],label=read_len$num_reads[m], vjust=0)
  
  p_sequence_length <- p1 +
    ggplot2::scale_x_continuous(breaks= seq(0,nrow(read_len),by=10),
                                limits = c(l,NA)) +
    ggplot2::labs(x = "Sequence length",
                  y = "Number of reads with sequence length",
                  title = "Sequence length distribution")
  
  p_sequence_length_log <- p1 +
    ggplot2::scale_y_continuous(trans='log10')+
    ggplot2::scale_x_continuous(breaks= seq(0,nrow(read_len),by=10),
                                limits = c(l,NA)) +
    ggplot2::labs(x = "Sequence length",
                  y = "Number of reads with sequence length",
                  title = "Logarithmic axis to show all reads")
  
  tx<-ggpubr::ggarrange(p_sequence_length,p_sequence_length_log,ncol=1,nrow=2)
  if(!is.na(output_file)){ggplot2::ggsave(file=output_file,tx,
                                          device='pdf',
                                          width=11.7, height=11.7,units='in')}
  return(tx)
}

read_lengths<-function(files=NULL, output=''){
  # Detailed plot of read length
  lengths <- matrix(ncol=2,nrow = 0)
  
  for (i in files) {
    # Fetch lengths
    fastq <- suppressMessages(seqTools::fastqq(i))
    read_len <- seqTools::seqLenCount(fastq)
    # Add rows so master matrix is right length
    if (nrow(lengths)<nrow(read_len)) {
      difference <- nrow(read_len)-nrow(lengths)
      add_rows   <- cbind((1:difference)+nrow(lengths),
                          rep(0,difference))
      lengths<-rbind(lengths,
                     add_rows)
    }
    # Tally up read numbers
    lengths[,2]<-lengths[,2]+read_len[,1]
  }
  
  # All files done? Plot:
  lengths <- as.data.frame(lengths)
  colnames(lengths)<-c('read_length','num_reads')
  plot_read_length(lengths, output_file = output)
}


##### Define main #####
main<-function(opt=NULL){
  #### Setup ####
  
  # Check all requirements
  Install_Check()
  
  # Check the arguments (with basic error conditions)
  if(is.null(opt$fastq_folder) | 
     is.null(opt$fw_primer) |
     is.null(opt$rv_primer)){
    warning('Missing required file. Please check given input')
    print(opt$fastq_folder)
    print(opt$fw_primer)
    print(opt$rv_primer)
    input_files()
  } else if(!dir.exists(opt$fastq_folder)){
    warning('Given folder cannot be found. Please check file path')
    fastqfolder()
  }
  
  # Pull all fastq files in the folder
  setwd(opt$fastq_folder)
  files <- list.files()
  
  #### Generate report folder ####
  if (!dir.exists('./Reports')) {
    dir.create('./Reports')}
  
  #### Data is paired end #####
  if(!is.null(opt$rv_primer)){
    
    #### Fetch fastq files ####
    # Put all forward reads in one object and the reverse reads in another
    fw_reads <- c(files[grep('.*\\_1.fq',files)],files[grep('.*\\_1.fastq',files)])
    rv_reads <- c(files[grep('.*\\_2.fq',files)],files[grep('.*\\_2.fastq',files)])
    
    #### Filter N reads ####
    # Make a new directory for N filtering
    if (!dir.exists('./Removed_N_Reads')) {
      dir.create('./Removed_N_Reads')}
    dir_fw_N_removed <- file.path('./Removed_N_Reads', basename(fw_reads))
    dir_rv_N_removed <- file.path('./Removed_N_Reads', basename(rv_reads))
    
    # Files may be overwritten
    if (opt$overwrite==T) {
      # Remove reads with one or more N
      message('Filtering out all reads with N bases\n')
      removed.N<-dada2::filterAndTrim(fwd  = fw_reads,
                                      filt = dir_fw_N_removed,
                                      rev  = rv_reads,
                                      filt.rev = dir_rv_N_removed,
                                      maxN = 0, 
                                      multithread = TRUE,
                                      verbose = T)
      # Save results
      percentages <- round(100*removed.N[,2]/removed.N[,1],2)
      report <- data.frame(cbind('Unfiltered_Total' = removed.N[,1],
                                 'After_N_Removal' = paste(removed.N[,2],
                                                           ' (',percentages,'%)',sep = '')))
      
      write.table(report,'./Reports/Report_Filtered_Reads.txt',sep='\t',quote=F)
      
      # Print overview to the terminal
      message('Removal of reads with one N or more resulted in:')
      print(cbind('Lowest_%_Remaining' = min(percentages),
                  'Highest_%_Remaining'= max(percentages),
                  'Median_%_Remaining' = round(median(percentages),2),
                  'Mean_%_Remaining'   = round(mean(percentages),2)))
      
      # Check if we lost any samples
      filtered_files <- list.files('./Removed_N_Reads')
      dir_fw_N_removed <- c(filtered_files[grep('.*\\_1.fq',filtered_files)],
                            filtered_files[grep('.*\\_1.fastq',filtered_files)])
      dir_fw_N_removed <- paste('./Removed_N_Reads/',dir_fw_N_removed,sep='')
      
      dir_rv_N_removed <- c(filtered_files[grep('.*\\_2.fq',filtered_files)],
                            filtered_files[grep('.*\\_2.fastq',filtered_files)])
      dir_rv_N_removed <- paste('./Removed_N_Reads/',dir_rv_N_removed,sep='')
      lost <- fw_reads[which(!fw_reads %in% basename(dir_fw_N_removed))]
      
      # And report losses if necessary
      if (length(lost)>0) {
        message('\nThe following samples have been filtered out completely:')
        print(lost)
      } else {
        message('\nNo samples were completely filtered out')
      }  
    }else if (!file.exists(dir_fw_N_removed[1])){
      warning('No N filtered folder available, primer filtering wont be done either')
    } else{message('N filtered files were not overwritten\n')}
    
    #### Only do primer work if necessary
    # Create a folder for primer removed files
    if (!dir.exists('./Removed_Primers')) {
      dir.create('./Removed_Primers')}
    dir_fw_primer_removed <- file.path('./Removed_Primers', basename(dir_fw_N_removed))
    dir_rv_primer_removed <- file.path('./Removed_Primers', basename(dir_rv_N_removed))
    
    if ((opt$overwrite==T | !file.exists(dir_fw_primer_removed[1])) & file.exists(dir_fw_N_removed[1])) {
      
      #### Find given primers ####
      # Get all possible primer variants
      fw_options <- allOrients(opt$fw_primer)
      rv_options <- allOrients(opt$rv_primer)
      
      # Check which orientation is correct
      message('Primer matches found:')
      hits <-
        rbind(FWD.ForwardReads = sapply(fw_options, primerHits, fn = dir_fw_N_removed[[1]]), 
              FWD.ReverseReads = sapply(fw_options, primerHits, fn = dir_rv_N_removed[[1]]), 
              REV.ForwardReads = sapply(rv_options, primerHits, fn = dir_fw_N_removed[[1]]), 
              REV.ReverseReads = sapply(rv_options, primerHits, fn = dir_rv_N_removed[[1]]))
      print(hits)
      message('\n')
      
      # Do we need to reorient any primers based on this search?
      fw_match <- c(grep(max(hits[1,]),hits[1,]),
                    grep(max(hits[2,]),hits[2,]))
      if(fw_match[1]!=1 & fw_match[2]!=4){
        message('Forward primer orientation was force set to ',
                names(fw_options[fw_match[1]]),'\n')
        opt$fw_primer<-fw_options[fw_match[1]]
      }
      
      rv_match <- c(grep(max(hits[3,]),hits[3,]),
                    grep(max(hits[4,]),hits[4,]))
      if(rv_match[1]!=4 & rv_match[2]!=1){
        message('Reverse primer orientation was force set to ',
                names(fw_options[rv_match[2]]),'\n')
        opt$fw_primer<-rv_options[rv_match[2]]
      }
      
      #### Remove found primers ####
      # Check if CUTADAPT is available
      present<-system2(file.path(opt$cutadapt), args = '--version',stdout = TRUE)
      
      condition <- length(present)>1
      if (condition) {
        print('Cutadapt module could not be reached!')
        cutadapt()
      }
      
      # Remove primers with CUTADAPT and full file paths
      Trim_fw <- paste("-g",opt$fw_primer,"-a",dada2::rc(opt$rv_primer))
      Trim_rv <- paste("-G",opt$rv_primer,"-A",dada2::rc(opt$fw_primer))
      base_path <- getwd()
      
      for (i in 1:length(dir_fw_primer_removed)) {
        message(basename(dir_fw_primer_removed[i]))
        message(basename(dir_rv_primer_removed[i]))
        system2(opt$cutadapt, args=c(Trim_fw,
                                     Trim_rv,
                                     "-n",2,
                                     "-o", sub('^\\.',base_path,dir_fw_primer_removed[i]),
                                     "-p", sub('^\\.',base_path,dir_rv_primer_removed[i]),
                                     sub('^\\.',base_path,dir_fw_N_removed[i]),
                                     sub('^\\.',base_path,dir_rv_N_removed[i]),
                                     "--minimum-length", 20,
                                     "--report=minimal",
                                     "--cores",0))
      }
      # Sanity check: no files were randomly lost 
      lost_fw <- which(!file.exists(dir_fw_primer_removed))
      lost_rv <- which(!file.exists(dir_rv_primer_removed))
      if (length(lost_fw)!=0) {
        warning(paste('Lost the following file:',dir_fw_primer_removed[lost_fw]))
      } else if (length(lost_rv!=0)) {
        warning(paste('Lost the following file:',dir_rv_primer_removed[lost_rv]))
      } else{message('\nNo files were lost')}
      
      # Sanity check: are the primers truly removed
      message('\n Check once more (Everything should be 0:')
      hits <-
        rbind(FWD.ForwardReads = sapply(fw_options, primerHits, fn = dir_fw_primer_removed[[1]]), 
              FWD.ReverseReads = sapply(fw_options, primerHits, fn = dir_rv_primer_removed[[1]]), 
              REV.ForwardReads = sapply(rv_options, primerHits, fn = dir_fw_primer_removed[[1]]), 
              REV.ReverseReads = sapply(rv_options, primerHits, fn = dir_rv_primer_removed[[1]]))
      print(hits)
      message('\n')
    } else if (!file.exists(dir_fw_N_removed[1])){
      warning('Primer filtering was not done due to overwrite=T. Set it to false in SLURM instead') 
    } else{message('Primer filtered files already exist and were not overwritten\n')}
    
    #### Plot quality ####
    # Now produce some quality plots for the user to help with parameter settings
    if (length(dir_fw_primer_removed)>12){
      files_subset<- sample(1:length(dir_fw_primer_removed),12)
    } else {files_subset <- 1:length(dir_fw_primer_removed)}
    base_path <- getwd()
    
    message('Plotting max 12 random forward read files\n')
    Qual_Plot_Name <- paste(getwd(),"/Reports/Forward_Read_Quality.pdf",sep="")

    tx<-dada2::plotQualityProfile(sub('^\\.',base_path,dir_fw_primer_removed[files_subset]))
    ggplot2::ggsave(Qual_Plot_Name,device='pdf',plot=tx,
                    width=11.7, height=11.7,units='in')
    
    message('Plotting max 12 random reverse read files\n')
    Qual_Plot_Name <- paste(getwd(),"/Reports/Reverse_Read_Quality.pdf",sep="")

    tx<-dada2::plotQualityProfile(sub('^\\.',base_path,dir_rv_primer_removed[files_subset]))
    ggplot2::ggsave(Qual_Plot_Name,device='pdf',plot=tx,
                    width=11.7, height=11.7,units='in')
    
    #### Plot read length ####
    # With previously set functions
    read_lengths(dir_fw_primer_removed,
                 './Reports/FW_Read_Length_Distribution.pdf')
    
    read_lengths(dir_rv_primer_removed,
                 './Reports/RV_Read_Length_Distribution.pdf')
    
    
    #### Benchmark filter thresholds ####
    # Test various filter thresholds to see how many reads are lost
    if (!dir.exists('./Parameter_Test')) {
      dir.create('./Parameter_Test')}
    dir_fw_par_test <- file.path('./Parameter_Test', basename(dir_fw_primer_removed))
    dir_rv_par_test <- file.path('./Parameter_Test', basename(dir_rv_primer_removed))
    
    message('Testing filter thresholds on the subset of files\n')
    maxEE.0.5<-dada2::filterAndTrim(fwd  = dir_fw_primer_removed[files_subset],
                                    filt = dir_fw_par_test[files_subset],
                                    rev  = dir_rv_primer_removed[files_subset],
                                    filt.rev = dir_rv_par_test[files_subset],
                                    maxEE = 0.5, 
                                    multithread = TRUE)
    
    maxEE.1<-dada2::filterAndTrim(fwd  = dir_fw_primer_removed[files_subset],
                                  filt = dir_fw_par_test[files_subset],
                                  rev  = dir_rv_primer_removed[files_subset],
                                  filt.rev = dir_rv_par_test[files_subset],
                                  maxEE = 1, 
                                  multithread = TRUE)
    
    maxEE.2<-dada2::filterAndTrim(fwd  = dir_fw_primer_removed[files_subset],
                                  filt = dir_fw_par_test[files_subset],
                                  rev  = dir_rv_primer_removed[files_subset],
                                  filt.rev = dir_rv_par_test[files_subset],
                                  maxEE = 2, 
                                  multithread = TRUE)
    
    maxEE.3<-dada2::filterAndTrim(fwd  = dir_fw_primer_removed[files_subset],
                                  filt = dir_fw_par_test[files_subset],
                                  rev  = dir_rv_primer_removed[files_subset],
                                  filt.rev = dir_rv_par_test[files_subset],
                                  maxEE = 3, 
                                  multithread = TRUE)
    
    phyx_on<-dada2::filterAndTrim(fwd  = dir_fw_primer_removed[files_subset],
                                  filt = dir_fw_par_test[files_subset],
                                  rev  = dir_rv_primer_removed[files_subset],
                                  filt.rev = dir_rv_par_test[files_subset],
                                  rm.phix = T, 
                                  multithread = TRUE)
    
    trunq <- dada2::filterAndTrim(fwd  = dir_fw_primer_removed[files_subset],
                                  filt = dir_fw_par_test[files_subset],
                                  rev  = dir_rv_primer_removed[files_subset],
                                  filt.rev = dir_rv_par_test[files_subset],
                                  minLen = 20,
                                  truncQ = 11, 
                                  multithread = TRUE)
    
    
    #### Save reports ####
    filter_effects <- cbind('Reads.In' = maxEE.0.5[,1],
                            'maxEE = 0.5' = paste(maxEE.0.5[,2], sep = '', ' (',
                                                  round(100*maxEE.0.5[,2]/maxEE.0.5[,1]),'%)'),
                            'MaxEE = 1' = paste(maxEE.1[,2], sep = '', ' (',
                                                round(100*maxEE.1[,2]/maxEE.0.5[,1]),'%)'),
                            'MaxEE = 2' = paste(maxEE.2[,2], sep = '', ' (',
                                                round(100*maxEE.2[,2]/maxEE.0.5[,1]),'%)'),
                            'MaxEE = 3' = paste(maxEE.3[,2], sep = '', ' (',
                                                round(100*maxEE.3[,2]/maxEE.0.5[,1]),'%)'),
                            'phiX genome removal' = paste(phyx_on[,2], sep = '', ' (',
                                                          round(100*phyx_on[,2]/maxEE.0.5[,1]),'%)'))
    message('\nNumber of reads remaining after the following filter thresholds:')
    print(filter_effects)
    
    message('\nThis table has also been saved as Filter_Tests.txt')
    write.table(filter_effects, file = './Reports/Filter_Tests.txt',sep='\t',quote=F)
    
    #### Remove N filtered files: they are not used downstream ####
    unlink('./Parameter_Test',recursive=T)
    if (opt$delete) {
      message('Removing excess directories to save space')
      unlink('./Removed_N_Reads',recursive=T)
    }
    
    # Data is not paired end
  } else{
    warning('Given data is not paired end!\nClosing pipeline.')
  }
}




