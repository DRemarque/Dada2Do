#### Collect all functions used in main ####
Argument_Parsing<-function(show_help=F){
  # Index Input
  option_list = list(
    optparse::make_option(c('--fastq_folder'), type='character', metavar = 'character', default=NULL,
                          help='The directory with all raw fastq files and filtered files in subfolders'),
    optparse::make_option(c('--novaseq'), type='logical', metavar = 'logical', default=T,
                          help='Whether or not it was sequenced with NovaSeq'),
    
    optparse::make_option(c('--otu_similarity'), type='numeric', metavar = 'numeric', default=0,
                          help='0.97 for 97% similarity clusters. Set to 0 to receive no OTU and treat all sequences as unique'),
    optparse::make_option(c('--otu_abundance'), type='logical', metavar = 'logical', default=T,
                          help='Is the otu representative sequence the most abundant sequence (T)\n
                          or is it the consensus of all sequences in that otu (F)'),
    

    optparse::make_option(c('--maxEE'), type='numeric', default=1,metavar = 'numerical',
                          help='Maximal estimated error of a read'),
    optparse::make_option(c('--truncLen_fw'), type='integer', metavar = 'integer', default=0,
                          help='Forward read length after which to truncate'),
    optparse::make_option(c('--truncLen_rv'), type='integer', metavar = 'integer', default=0,
                          help='Reverse read length after which to truncate'),
    
    optparse::make_option(c('--minLen'), type='integer', metavar = 'integer', default=20,
                          help='Minimum length of a read after truncating'),
    optparse::make_option(c('--truncQ'), type='integer', metavar = 'integer', default=2,
                          help='Minimal basequality before truncating'),
    
    optparse::make_option(c('--pool_samples'), type='character', metavar = 'character', default=FALSE,
                          help='Whether to pool samples: TRUE, FALSE or PSEUDO. The last is an intermediate option'),
    
    optparse::make_option(c('--database'), type='character', metavar = 'character', default='./dummy.fa',
                          help='Database file name'),
    
    optparse::make_option(c('--overwrite'), type='logical', metavar = 'logical', default=T,
                          help='Wheter previously made files should be overwritten')
  )
  opt_parser <- optparse::OptionParser(option_list=option_list)
  if(show_help){
    optparse::print_help(opt_parser)
  } else{
    opt <- optparse::parse_args(opt_parser)
    return(opt)  
  }
}#Argument_Parsing

Install_Check<-function(){
  installed <- installed.packages()
  cran <- 'https://mirror.lyrahosting.com/CRAN/'
  if (!'BiocManager' %in% installed) {
    message('Bioconductor not found. Installing now\n')
    suppressMessages(install.packages('BiocManager',repos=cran))
  }
  if (!'dada2' %in% installed) {
    message('dada2 not found. Installing now\n')
    suppressMessages(BiocManager::install('dada2'))
  }
  if (!'optparse' %in% installed) {
    message('optparse not found. Installing now\n')
    suppressMessages(install.packages('optparse',repos=cran))
  }
  if (!'dplyr' %in% installed) {
    message('optparse not found. Installing now\n')
    suppressMessages(install.packages('optparse',repos=cran))
  }
  if (!'DECIPHER' %in% installed) {
    message('DECIPHER not found in correct version. Installing now\n')
    suppressMessages(BiocManager::install('DECIPHER'))
  } 
  if (!'Biostrings' %in% installed) {
    message('Biostrings not found. Installing now\n')
    suppressMessages(BiocManager::install('Biostrings'))
  }
  if (!'msa' %in% installed) {
    message('MSA not found. Installing now\n')
    suppressMessages(BiocManager::install('msa'))
  }
  message('All dependencies available\n')
  
}#Install_Check

log_filter_results <- function(general_file = NULL,
                               pool_specific_file = NULL,
                               new_report = NULL,
                               col_ID = NULL){
  
  # We are writing in the general file
  if(is.null(pool_specific_file)){
    if(!file.exists(general_file)){
      # warn user
      warning('The filter report of the previous pipeline step was not found. \n Creating a new one...')
      if (!dir.exists('./Reports')) {
        dir.create('./Reports')
      }
      # Write and save
      write.table(new_report, file = general_file,quote=F,sep='\t')
      
    } else{
      # read table
      report <- read.table(general_file, check.names = F,sep='\t')
      
      # Overwrite previous filter and trim column?
      if(TRUE %in% grepl(col_ID,colnames(report))){
        ind <- grep(col_ID,colnames(report))
        report <- report[,-ind]} 
      
      # Combine and save
      report_merged <- merge(report, new_report[,-1], by=0,all=T)
      report_merged<- report_merged[,-1]
      colnames(report_merged)<-c(colnames(report),
                                 colnames(new_report)[-1])
      rownames(report_merged)<-rownames(report)
      write.table(report_merged, file = general_file,quote=F,sep='\t')
      
    } # overwrite report table?
    
    # We are writing in the pool specific file 
  } else{
    # OPTION 1: No reports exist at all
    if(!file.exists(general_file) &
       !file.exists(pool_specific_file)){
      # warn user of missing report file
      warning('The filter report of the previous pipeline steps was not found. 
              \n Creating a new one...')
      if (!dir.exists('./Reports')) {
        dir.create('./Reports')
      }
      # Write and save
      write.table(new_report, file = pool_specific_file,quote=F,sep='\t')
      
      # OPTION 2: General file exists to hitch onto
    } else if(file.exists(general_file) &
              !file.exists(pool_specific_file)){
      # read table
      report <- read.table(general_file, check.names = F,sep='\t')
      
      # Combine and save
      report <- merge(report,new_report, by=0,all=T)
      write.table(report, file = pool_specific_file,quote=F,sep='\t')
      
      # OPTION 3: Last resort, general file does not exist but pool specific one does
    } else{
      # read table
      report <- read.table(pool_specific_file, check.names = F,sep='\t')
      
      # Overwrite previous filter and trim column?
      if(TRUE %in% grepl(col_ID,colnames(report))){
        report[,grep(col_ID,colnames(report))] <- 
          report[,-grep(col_ID,colnames(report))]} 
      
      # Combine and save
      report <- merge(report,new_report, by=0,all=T)
      write.table(report, file = pool_specific_file,quote=F,sep='\t')
    } # overwrite report table?
    
  }
}#log_filter_results

loessErrfun_mod4 <- function(trans) {
  # This is a function that LearnError can use to set Loess
  # It sets weights, span and degree while also enforcing monotonicity
  
  # Ensure correct input and create matrix to be filled
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  
  # Do all comparisons between nucleotide transversions (i.e. A>C and C>A)
  for(nti in c('A','C','G','T')) {
    for(ntj in c('A','C','G','T')) {
      # A comparison is not necessary when A>A or T>T: there is no change
      if(nti != ntj) {
        # Title of the comparison (i.e. A2C or G2T)
        errs <- trans[paste0(nti,'2',ntj),]
        # Total transition rate of one nucleotide to all others
        tot <- colSums(trans[paste0(nti,'2',c('A','C','G','T')),])
        # 1 psuedo-count for each err, but if tot=0 will give NA
        rlogp <- log10((errs+1)/tot) 
        rlogp[is.infinite(rlogp)] <- NA
        # Put all of the above in one frame
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        
        # original
        # ###! mod.lo <- loess(rlogp ~ q, df, weights=errs) ###!
        # mod.lo <- loess(rlogp ~ q, df, weights=tot) ###!
        # #        mod.lo <- loess(rlogp ~ q, df)
        
        # jonalim's solution
        # https://github.com/benjjneb/dada2/issues/938
        
        # Perform the customised Loess model
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),degree = 1, span = 0.95)
        
        # Use the model to predict the error rates
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        
        # Store the predictions in the overarching matrix
        est <- rbind(est, 10^pred)
      } # if(nti != ntj)
    } # for(ntj in c('A','C','G','T'))
  } # for(nti in c('A','C','G','T'))
  
  # HACKY (enforces fixed max and min rates)
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  # Enforce monotonicity
  # https://github.com/benjjneb/dada2/issues/791
  estorig <- est
  est <- est %>%
    data.frame() %>%
    dplyr::mutate_all(dplyr::funs(dplyr::case_when(. < X40 ~ X40,
                                                   . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  # Expand the error matrix with the self-transition probabilities
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c('A','C','G','T'), each=4), '2', c('A','C','G','T'))
  colnames(err) <- colnames(trans)
  # Return
  return(err)
}#loessErrfun_mod4

converge_OK<-function(error_model){
  # Were the rates close to convergence?
  conv <- dada2:::checkConvergence(error_model)
  conv <- conv[(length(conv)-3):length(conv)]
  # No 
  if (max(dist(conv))>0.1){
    message('Error rate differences: ',round(max(dist(conv))),2)
    message('Non convergence in error learning!\n',
             'Something seems wrong with the error rates of your data\n',
             'Please check the quality and use pipeline results with apprehension!!!')
  } else{
    message('Error learning reached decent convergence\n',
            'Still, do check if the error rate plots show good fits')
  }# rates ok?
}#converge_OK

getN <- function(x) sum(dada2::getUniques(x))#getN

# DECIPHER's old IdCluster function (no longer available :( )
source('./DECIPHER_OLD_IdClusters.R')

##### Define main #####
main<-function(opt=NULL){
  #### Housekeeping ####
  # Check all requirements
  Install_Check()
  
  # Check the arguments (with basic error conditions)
  if(is.null(opt$fastq_folder)){
    warning('Missing required file. Please check given input')
    print(opt$fastq_folder)
    input_files()
  } else if(!dir.exists(opt$fastq_folder)){
    warning('Given folder cannot be found. Please check file path')
    print(opt$fastq_folder)
    fastqfolder()
  } else if(!dir.exists(paste(sub('/$','',opt$fastq_folder),
                              '/Removed_Primers',sep=''))){
    warning('No primer free fastq folder can be found. Please check file path or perform the previous pipeline step first.\n',
            paste(sub('/$','',opt$fastq_folder),'/Removed_Primers',sep=''))
    fastqNfolder()
  } else if(opt$pool_samples!=TRUE   &
            opt$pool_samples!='TRUE' &
            opt$pool_samples!=FALSE  &
            opt$pool_samples!='FALSE'&
            opt$pool_samples!='PSEUDO'){
    warning('Invalid pool argument: use TRUE, FALSE or PSEUDO')
    pool_argument()
  } # input check
  
  # Pull all fastq files in the folder
  setwd(opt$fastq_folder)
  files <- list.files('./Removed_Primers')
  
  # Load pipe symbol
  suppressMessages(library('dplyr'))
  
  #### Fetch fastq files ####
  # Put all forward reads in one object and the reverse reads in another
  dir_fw_primer_removed <- c(files[grep('.*\\_1.fq',files)],files[grep('.*\\_1.fastq',files)])
  dir_fw_primer_removed <- paste('./Removed_Primers/',dir_fw_primer_removed,sep='')
  dir_rv_primer_removed <- c(files[grep('.*\\_2.fq',files)],files[grep('.*\\_2.fastq',files)])
  dir_rv_primer_removed <- paste('./Removed_Primers/',dir_rv_primer_removed,sep='')
  
  #### Trimming and filtering ####
  # Setup a new folder
  if (!dir.exists('./Filtered_Reads')) {
    dir.create('./Filtered_Reads')}
  
  dir_fw_filtered <- file.path('./Filtered_Reads', basename(dir_fw_primer_removed))
  dir_rv_filtered <- file.path('./Filtered_Reads', basename(dir_rv_primer_removed))
  
  # Trim only if files may be overwritten
  if (opt$overwrite==T | !file.exists(dir_fw_filtered[1])) {
    message('Filtering and trimming...\n')
    filtered_report <- dada2::filterAndTrim(fwd  = dir_fw_primer_removed,
                                            filt = dir_fw_filtered,
                                            rev  = dir_rv_primer_removed,
                                            filt.rev = dir_rv_filtered,
                                            truncLen = c(opt$truncLen_fw,
                                                         opt$truncLen_rv),
                                            maxEE  = opt$maxEE,
                                            minLen = opt$minLen,
                                            truncQ = opt$truncQ,
                                            maxN   = 0, 
                                            verbose = T,
                                            multithread = TRUE)
    # Collect results
    percentages <- round(100*filtered_report[,2]/filtered_report[,1],2)
    colname <- paste('truncQ:',opt$truncQ,
                     ' minLen:',opt$minLen,
                     ' MaxEE:',opt$maxEE, sep = '')
    filtered_report[,2] <- paste(filtered_report[,2], ' (', percentages,'%)',sep='')
    colnames(filtered_report)[2] <-colname
    
    # Log to correct file
    log_filter_results(general_file = './Reports/Report_Filtered_Reads.txt',
                       new_report = filtered_report,
                       col_ID = 'MaxEE')
    
    # Print general results to the output log
    message('Read filtering and trimming resulted in:')
    print(cbind('Lowest_%_Remaining' = min(percentages),
                'Highest_%_Remaining'= max(percentages),
                'Median_%_Remaining' = round(median(percentages),2),
                'Mean_%_Remaining'   = round(mean(percentages),2)))
    
    # Check if we lost any samples
    filtered_files <- list.files('./Filtered_Reads')
    dir_fw_filtered <- c(filtered_files[grep('.*\\_1.fq',filtered_files)],
                         filtered_files[grep('.*\\_1.fastq',filtered_files)])
    dir_fw_filtered <- paste('./Filtered_Reads/',dir_fw_filtered,sep='')
    dir_rv_filtered <- c(filtered_files[grep('.*\\_2.fq',filtered_files)],
                         filtered_files[grep('.*\\_2.fastq',filtered_files)])
    dir_rv_filtered <- paste('./Filtered_Reads/',dir_rv_filtered,sep='')
    lost_fw <- dir_fw_primer_removed[which(!basename(dir_fw_primer_removed) %in% basename(dir_fw_filtered))]
    lost_rv <- dir_rv_primer_removed[which(!basename(dir_rv_primer_removed) %in% basename(dir_rv_filtered))]
    
    # And report losses if necessary (sanity check)
    if (length(lost_fw)>0) {
      message('\nThe following forward samples have been lost:')
      print(lost_fw)
    } else if(length(lost_rv)>0){
      message('\nThe following reverse samples have been lost:')
      print(lost_rv)
    } else {
      message('\nNo samples were completely filtered out')
    }
  }# overwrite | !file.exists filtered reads
  
  #### Creating updated quality plots ####
  # Now produce some quality plots for control
  if (length(dir_fw_filtered)>12){
    files_subset<- sample(1:length(dir_fw_filtered),12)
  } else {files_subset <- 1:length(dir_fw_filtered)} #subset samples
  
  Qual_Plot_Name <- paste(getwd(),'/Reports/Filtered_Forward_Read_Quality.pdf',sep='')
  
  if (opt$overwrite == T | !file.exists(Qual_Plot_Name)) {
    message('Plotting max 12 random forward read files\n')
    
    tx<-dada2::plotQualityProfile(dir_fw_filtered[files_subset])
    ggplot2::ggsave(Qual_Plot_Name,device='pdf',plot=tx,
                    width=11.7, height=11.7,units='in')
  }# overwrite fw qual plot
  
  Qual_Plot_Name <- paste(getwd(),'/Reports/Filtered_Reverse_Read_Quality.pdf',sep='')
  
  if (opt$overwrite == T | !file.exists(Qual_Plot_Name)) {
    message('Plotting max 12 random reverse read files\n')
    
    tx<-dada2::plotQualityProfile(dir_rv_filtered[files_subset])
    ggplot2::ggsave(Qual_Plot_Name,device='pdf',plot=tx,
                    width=11.7, height=11.7,units='in')
  } # overwrite rv qual plot
  
  #### Make error models if not yet made ####
  if (opt$overwrite | 
      !file.exists('./Reports/Forward_Error_Model.rds') |
      !file.exists('./Reports/Reverse_Error_Model.rds'))  {

    #### Error learning for MiSeq/HiSeq ####
    if(!opt$novaseq){
      message('\nError learning for MiSeq/HiSeq...\nset.seed=21\nUsing ',
              RcppParallel::defaultNumThreads(),' cores')
      set.seed(21)
      # Forward read learning
      message('\nForward reads: found ',length(dir_fw_filtered),' samples\n')
      fw_error <- dada2::learnErrors(dir_fw_filtered,
                                     nbases=1e8,
                                     MAX_CONSIST = 15,
                                     multithread=TRUE, randomize = TRUE)
      saveRDS(fw_error,'./Reports/Forward_Error_Model.rds')
      # Non convergence after 15 steps?
      converge_OK(fw_error)
      
      # Reverse read learning
      message('\nReverse reads: found ',length(dir_rv_filtered),' samples\n')
      rv_error <- dada2::learnErrors(dir_rv_filtered,
                                     nbases=1e8,
                                     MAX_CONSIST = 15,
                                     multithread=TRUE, randomize = TRUE)
      saveRDS(rv_error,'./Reports/Reverse_Error_Model.rds')
      # Non convergence after 15 steps?
      converge_OK(rv_error)
      
      #### Error learning for NovaSeq ####
    }else{
      message('\nError learning for Novaseq...\nset.seed=21\nUsing ',
              RcppParallel::defaultNumThreads(),' cores')
      set.seed(21)
      # Forward read learning
      message('Forward reads: found ',length(dir_fw_filtered),' samples\n')
      fw_error <- dada2::learnErrors(dir_fw_filtered,
                                     nbases = 1e8,
                                     errorEstimationFunction = loessErrfun_mod4, 
                                     randomize = T,
                                     MAX_CONSIST = 15,
                                     multithread =T,
                                     verbose = TRUE)
      saveRDS(fw_error,'./Reports/Forward_Error_Model.rds')
      # Non convergence after 15 steps?
      converge_OK(fw_error)
      
      # Reverse read learning
      message('Reverse reads: found ',length(dir_rv_filtered),' samples\n')
      rv_error <- dada2::learnErrors(dir_rv_filtered,
                                     nbases = 1e8,
                                     errorEstimationFunction = loessErrfun_mod4, 
                                     randomize = T,
                                     MAX_CONSIST = 15,
                                     multithread =T,
                                     verbose = TRUE)
      saveRDS(rv_error,'./Reports/Reverse_Error_Model.rds')
      # Non convergence after 15 steps?
      converge_OK(rv_error)
      
    } # novaseq | miseq learnError
    
    # Save the error models so they don't have to be remade
    message('Error models made and saved in ./Reports folder')
  } else{
    message('Old error models found in ./Reports folder. Using those!\n')
    fw_error <- readRDS('./Reports/Forward_Error_Model.rds')
    rv_error <- readRDS('./Reports/Reverse_Error_Model.rds')
  } # if overwrite | !file.exists 
  
  #### Show the error models ####
  if(opt$overwrite | 
     !file.exists('./Reports/Forward_Error_Model.pdf') |
     !file.exists('./Reports/Reverse_Error_Model.pdf')){
    message('Plotting error models\n')
    
    tx<-dada2::plotErrors(fw_error, nominalQ=T)
    ggplot2::ggsave('./Reports/Forward_Error_Model.pdf',
                    device='pdf',plot=tx,
                    width=11.7, height=11.7,units='in')
    
    tx<-dada2::plotErrors(rv_error, nominalQ=T)
    ggplot2::ggsave('./Reports/Reverse_Error_Model.pdf',
                    device='pdf',plot=tx,
                    width=11.7, height=11.7,units='in')

  } # overwrite plots
  
  #### Dereplicate and merge reads ####
  # Only if necessary!
  raw_sequence_file <- paste('./Reports/Pool_',opt$pool_samples,
                             '_Raw_Sequence_Table.rds',sep='')
  if(opt$overwrite |
     !file.exists(raw_sequence_file)){
    
    # How many samples?
    samples <- basename(dir_fw_filtered)
    
    # Pooled samples
    if(opt$pool_samples=='TRUE' | 
       opt$pool_samples==TRUE |
       opt$pool_samples=='PSEUDO'){
      message('Dereplicating and merging samples with pool = ',opt$pool_samples)
      pool_type <- ifelse((opt$pool_samples=='TRUE' | opt$pool_samples==TRUE) ,TRUE,'PSEUDO')
      
      dereplicated_fw <- dada2::derepFastq(dir_fw_filtered)
      dereplicated_rv <- dada2::derepFastq(dir_rv_filtered)
      
      names(dereplicated_fw) <- samples
      names(dereplicated_rv) <- samples
      
      error_corrected_fw <- dada2::dada(dereplicated_fw, err = fw_error, 
                                        multithread = T,
                                        pool = pool_type)
      
      
      error_corrected_rv <- dada2::dada(dereplicated_rv, err = rv_error, 
                                        multithread = T,
                                        pool = pool_type)
      
      names(error_corrected_fw) <- samples
      names(error_corrected_rv) <- samples
      
      merged_reads <- dada2::mergePairs(error_corrected_fw, dereplicated_fw,
                                        error_corrected_rv, dereplicated_rv)
      
    }else{
      message('Dereplicating and merging samples with pool = FALSE')
      merged_reads <- vector('list', length(dir_fw_filtered))
      error_corrected_fw <- vector("list", length(dir_fw_filtered))
      error_corrected_rv <- vector("list", length(dir_rv_filtered))
      
      names(merged_reads) <- samples
      names(error_corrected_fw) <- samples
      names(error_corrected_rv) <- samples
      
      for (sampleNR in 1:length(dir_fw_filtered)) {
        message('Processing ',sampleNR,' out of ',length(dir_fw_filtered),': ', names(merged_reads[sampleNR]))
        # Dereplicate
        dereplicated_fw <- dada2::derepFastq(dir_fw_filtered[sampleNR])
        dereplicated_rv <- dada2::derepFastq(dir_rv_filtered[sampleNR])
        
        # Error correct
        dada_fw <- dada2::dada(dereplicated_fw, err = fw_error, 
                               multithread = T,
                               pool = F)
        error_corrected_fw[[sampleNR]] <- dada_fw
        
        dada_rv <- dada2::dada(dereplicated_rv, err = rv_error, 
                               multithread = T,
                               pool = F)
        error_corrected_rv[[sampleNR]] <- dada_rv
        
        # Merge everything
        temp_merged <- dada2::mergePairs(dada_fw, dereplicated_fw,
                                         dada_rv, dereplicated_rv)
        
        merged_reads[[sampleNR]] <- temp_merged
      } # for loop
    } # pool type
    # Save just in case
    saveRDS(error_corrected_fw,paste0('./Reports/Pool_',opt$pool_samples,
                                      '_Error_Corrected_FW.rds'))
    saveRDS(error_corrected_rv,paste0('./Reports/Pool_',opt$pool_samples,
                                      '_Error_Corrected_RV.rds'))
    
    message('Done with denoising and merging\n')
    
    #### Raw ASV Table ####
    # Take the merged reads into a sequence table
    Sequence_Table_Raw <- dada2::makeSequenceTable(merged_reads)
    saveRDS(Sequence_Table_Raw, file = raw_sequence_file)
    
    #### Collect denoising and merging results ####
    ind <- which(basename(dir_fw_filtered) %in% samples)
    Denoised_Results <- data.frame(cbind('Denoised_fw' = 
                                           paste(sapply(error_corrected_fw[samples], getN),
                                                 ' (',
                                                 round(100*
                                                         sapply(error_corrected_fw[samples], getN)/
                                                         sapply(dir_fw_filtered[ind], getN),2), 
                                                 '%)',sep = ''),
                                         'Denoised_rv' = 
                                           paste(sapply(error_corrected_rv[samples], getN),
                                                 ' (',
                                                 round(100*
                                                         sapply(error_corrected_rv[samples], getN)/
                                                         sapply(dir_rv_filtered[ind], getN),2), 
                                                 '%)',sep = ''),
                                         'After_Merge' =
                                           sapply(merged_reads[samples], getN)))
    
    # Log denoising results
    log_filter_results(general_file = './Reports/Report_Filtered_Reads.txt',
                       pool_specific_file = paste0('./Reports/Pool_',
                                                  opt$pool_samples,
                                                  '_Report_Filtered_Reads.txt'),
                       new_report = Denoised_Results,
                       col_ID = 'Denoised_fw')
    
  }else{
    message('Old raw sequence table found, using this one\n',
            raw_sequence_file,'\n')
    Sequence_Table_Raw <- readRDS(raw_sequence_file)
  } # overwrite raw sequence table
  
  #### Chimera removal ####
  chimera_file <- paste('./Reports/Pool_',opt$pool_samples,
                        '_NonChim_Sequence_Table.rds',sep='')
  # Only if necessary
  if(opt$overwrite |
     !file.exists(chimera_file)){
    message('Removing chimera sequences...')
    # Remove
    Sequence_Table <- dada2::removeBimeraDenovo(Sequence_Table_Raw, 
                                                method = 'consensus',
                                                multithread = T)
    # Save chimera table
    saveRDS(Sequence_Table,chimera_file)
    write.csv(Sequence_Table,sub('rds$','csv',chimera_file))
    
    # Collect results
    nonChimResults <- data.frame('Non_Chimeric'=paste(rowSums(Sequence_Table),' (',
                                                      round(100*rowSums(Sequence_Table)/rowSums(Sequence_Table_Raw),2),
                                                      '%)',sep = ''))
    rownames(nonChimResults) <- rownames(Sequence_Table)
    
    # Save report
    log_filter_results(general_file = './Reports/Report_Filtered_Reads.txt',
                       pool_specific_file = paste('./Reports/Pool_',
                                                  opt$pool_samples,
                                                  '_Report_Filtered_Reads.txt',
                                                  sep=''),
                       new_report = nonChimResults,
                       col_ID = 'Non_Chimeric')
    
    # Show results
    message('The following percentage of total reads was non chimeric: ',
            round(100*sum(Sequence_Table)/sum(Sequence_Table_Raw),2),'%')
  }else{
    message('Old clean sequence table found, using this one\n',
            chimera_file,'\n')
    Sequence_Table <- readRDS(chimera_file)
  }# overwrite nonchim
  
  #### Excel save of filter report ####
  report_file <- paste0('./Reports/Pool_',
                        opt$pool_samples,
                        '_Report_Filtered_Reads.txt')
  if(file.exists(report_file)){
    report <- read.table(report_file, check.names = F,sep='\t')
    Sample_column <- grep('fq | fastq',report[1,])
    Sample_column <- ifelse
    
    Sample_Name <- gsub('\\_.\\.fq.*','',rownames(report))
    Sample_Name <- gsub('\\_.\\.fastq.*','',Sample_Name)
    report <- cbind('Sample_Name'=Sample_Name,
                    report)
    openxlsx::write.xlsx(report,file = paste0('./Reports/Pool_',
                                              opt$pool_samples,
                                              '_Final_Filter_Report.xlsx'))
  }
  
  #### Fungi only: OTU creation ####
  otu_file <- paste('./Reports/Pool_',opt$pool_samples,'_',
                    gsub('\\. | \\,','_',opt$otu_similarity),
                    '_OTU_Table.rds',sep='') 
  # Only if necessary
  if(opt$otu_similarity < 1 & 
     (opt$overwrite | !file.exists(otu_file))){
    message('Creating a new OTU table with ',opt$otu_similarity*100,'% clusters')
    message('\nUsing ',
            RcppParallel::defaultNumThreads(),' cores')
    # Collect info
    Sequence_Table <- readRDS(chimera_file)
    asv_sequences <- colnames(Sequence_Table)
    sample_names  <- rownames(Sequence_Table)
    dna <- Biostrings::DNAStringSet(asv_sequences)
    cutoff_value  <- 1-opt$otu_similarity # Usually 0.03
    
    ## Find clusters of ASVs to form the new OTUs
    #align    <- DECIPHER::AlignSeqs(dna, 
    #                                iterations = 0,
    #                                useStructures = F,
    #                                processors = RcppParallel::defaultNumThreads(),
    #                                verbose = T)
    alignment_file <- paste0('./Reports/Pool_',opt$pool_samples,
                             '_Alignment.rds')
    if(!file.exists(alignment_file)){
      align <- msa::msaMuscle(dna,type='dna',order='input',
                              maxiters=1,verbose=F,diags=T)
      align <- align@unmasked
      saveRDS(align, alignment_file)
      message('Saved sequence alignment')
    } else{
      align <- readRDS(alignment_file)
      message('Loaded previous ASV alignment')
    }
    
    message('\nGetting distance matrix')
    distance <- DECIPHER::DistanceMatrix(align, 
                                         processors = RcppParallel::defaultNumThreads(),
                                         verbose = T)
    message('\nIdentifying clusters')
    clusters <- IdClusters(distance, 
                           method = "complete",
                           cutoff = cutoff_value, 
                           processors = RcppParallel::defaultNumThreads(),
                           verbose = T)
    # Add a new column
    clusters <- data.frame(OTU = clusters[,1],
                           Unique_Sequence = asv_sequences)
    
    Merged_Sequence_Table <- Sequence_Table %>% 
      t %>%
      rowsum(clusters$OTU) %>%
      t
    
    # Fetch the most abundant sequence from each OTU for taxon ID
    OTU_sequences <- colnames(Merged_Sequence_Table)
    clusters$OTU_Representative <- rep(NA,nrow(clusters))
    clusters$Is_Representative  <- rep(NA,nrow(clusters))
    
    if(opt$otu_abundance==T){
      message('OTU representative sequence is most abundant sequence')
      for (otuNR in 1:length(OTU_sequences)) {
        # Get all sequences of this OTU
        index <- which(clusters$OTU == as.numeric(OTU_sequences[otuNR]))
        # Which is most abundant
        unique_sequences <- summary(as.factor(clusters$Unique_Sequence[index]))
        most_abundant <- names(which(unique_sequences==max(unique_sequences))[1])
        
        # Log this one
        OTU_sequences[otuNR] <- most_abundant
        
        # And also write it to the clusters object for backlogging
        clusters$OTU_Representative[index]<-most_abundant
        clusters$Is_Representative[index] <- sapply(clusters$Unique_Sequence[index],
                                                    function(x){
                                                      ifelse(x==most_abundant,TRUE,FALSE)
                                                    })
        
      }# for loop
    } else{
      message('OTU representative sequence is most consensus sequence')
      for (otuNR in 1:length(OTU_sequences)) {
        # Get all sequences of this OTU
        index <- which(clusters$OTU == as.numeric(OTU_sequences[otuNR]))
        sequence_set <- Biostrings::DNAStringSet(clusters$Unique_Sequence[index])
        # What is the consensus?
        OTU_consensus <- DECIPHER::ConsensusSequence(sequence_set)
        
        # Log this one 
        OTU_sequences[otuNR] <- as.character(OTU_consensus)
        
        # And also write it to the clusters object for backlogging
        clusters$OTU_Representative[index]<-as.character(OTU_consensus)
      } # for loop
    } # consensus or abundance
    # Save the OTU key table
    message('Saving sequence to OTU key and full OTU table in ./Reports')
    saveRDS(clusters,
            paste0('./Reports/','Pool_',opt$pool_samples,'_',
                   gsub('\\. | \\,','_',opt$otu_similarity),'_OTU_key.rds'))
    
    # Replace OTU names with most prominent sequence
    colnames(Merged_Sequence_Table) <- OTU_sequences
    
    # Save OTU table
    saveRDS(Merged_Sequence_Table, otu_file)
    write.csv(Merged_Sequence_Table,sub('rds$','csv',otu_file))
    
    # Report
    message('\nChanged ',ncol(Sequence_Table),' to ',ncol(Merged_Sequence_Table),' OTUs\n')
    
    # Overwrite ASV
    Sequence_Table <- Merged_Sequence_Table 
  } else if(opt$otu_similarity < 1){
    # Read and overwrite
    message('Old OTU table found, using this one\n',
            otu_file,'\n')
    Sequence_Table <- readRDS(otu_file)
  } # OTU?
  
  #### Assign taxonID ####
  # Get a basename
  database_name <- sub('\\.fa.*$','',basename(opt$database))
  database_name <- sub('\\.fasta.*$','',basename(database_name))
  
  # Generate a file name
  taxon_file <- paste('./Reports/Pool_',opt$pool_samples,
                      '_Taxon_Table_',database_name,
                      '.rds',sep='')
  
  # If the database exists and the file needs to be made
  if(file.exists(opt$database)){
    message('Using ', opt$database,': Which was last updated at ',
            file.info(opt$database)$ctime,'\n')
    
    message('Assigning taxons, this may take a while')
    # Assign taxons
    taxon_table <- dada2::assignTaxonomy(seqs = Sequence_Table,
                                         refFasta = opt$database,
                                         tryRC = T,
                                         multithread = T,
                                         verbose = T)
    # Clean up the table 
    
    
    # Save taxon
    message('Saving taxonomy data as ',sub('\\.rds$','',taxon_file))
    saveRDS(taxon_table, taxon_file)
    openxlsx::write.xlsx(data.frame(taxon_table), sub('rds$','xlsx',taxon_file),
                         colNames=TRUE, rowNames = TRUE, 
                         keepNA=TRUE, na.string='NA')

  }else{
    warning('Taxon database could not be found!\n
            Check the file path: ',opt$database,
            '\n and rerun the pipeline with overwrite=FALSE once the issue is fixed')
  }
}#main