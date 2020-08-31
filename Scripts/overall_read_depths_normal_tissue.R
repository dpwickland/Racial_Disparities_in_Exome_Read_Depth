#usage: Rscript --vanilla TCGA_overall_depths.R BRCA

#in output, each line corresponds to ONE individual

args<-commandArgs(TRUE)
.libPaths()
# test whether both arguments present
if (length(args)!=1) {
  stop("Usage example: Rscript --vanilla TCGA_overall_depths.R BRCA", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  CANCER_NAME <- toString(args[1])
}

#load libraries
.libPaths(c( .libPaths(), "/home/mayo/m187735/R", "/usr/local/biotools/rpackages/R-3.5.2-latest",  "/usr/local/biotools/r/R-3.5.2/lib64/R/library"))
library(plyr)
library(dplyr)
library(readr)
library(stringr)
setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)

#1. Load all clinical data ##CREATE THESE FILES BEFOREHAND
  all_clinical_data <- read.table(paste('/home/mayo/m187735/s212975.Wickland_Immunomics/TCGA_metadata/clinical/',CANCER_NAME,'_clinical_data.txt',sep=''),header=TRUE,stringsAsFactors = FALSE)
  
  #only need tumor stage and race
  all_clinical_data <- all_clinical_data[,c("submitter_id","tumor_stage","race")]


#2. Load immune data (from Thorsson et al, 2018)
  immune_data <- data.frame(read.table('/home/mayo/m187735/s212975.Wickland_Immunomics/TCGA_metadata/immune_data_Thorsson.txt', header=TRUE))
  #rename case ID column
  colnames(immune_data)[colnames(immune_data) == 'TCGA_Participant_Barcode'] <- 'submitter_id'
  #subset cancer type
  immune_data <- subset(immune_data, (TCGA_Study == CANCER_NAME))
  #remove superfluous column
  immune_data <- immune_data[-c(2)]

#5. Retrieve read depths of exome bams; read depths calculated by the script ***calc_read_depths_in_batches.txt***
  read_depths <- read.table(paste('/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/read_depths_F4_F1024_q20//',CANCER_NAME,'/all_read_depths.txt',sep=''),header=FALSE)
  colnames(read_depths) <- c("file_name","read_depth")
  read_depths$bam_full_path <- read_depths$file_name 
  read_depths$file_name <- sub(".*\\/", "\\1",read_depths$bam_full_path) #remove path from filename
  
  #Get metadata for all bams 
  exome_bams_metadata <- read.table(paste('/home/mayo/m187735/s212975.Wickland_Immunomics/TCGA_metadata/exome_bams/',CANCER_NAME,'_bams_metadata.txt',sep=''),header=TRUE)
  exome_bams_metadata$submitter_id <- sub("^([^-]*-[^-]*-[^-]*-[^-]*).*", "\\1",exome_bams_metadata$cases)
  
  #only interested in primary solid tumor
  exome_bams_metadata  <- subset(exome_bams_metadata, (tissue.definition=='Blood Derived Normal'))
  
  #extract tissue source site for exome
  #extract center from barcode; must remove everything before last hyphen
  exome_bams_metadata$center <- sub(".*\\-", "\\1",exome_bams_metadata$cases)
  
  #extract TSS code **from barcode**; remove everything after 2nd hyphen
  exome_bams_metadata$TSS_Code <- gsub("TCGA-","",sub("^([^-]*-[^-]*).*", "\\1",exome_bams_metadata$cases))
  
  #match TSS code with site ID (from https://gdc.cancer.gov/files/public/file/tcga_code_tables.zip) -- replace spaces, apostrophes with underscores first
  TSS <- read.table('/home/mayo/m187735/s212975.Wickland_Immunomics/TCGA_metadata/tissue_source_site.txt', header=TRUE)
  exome_bams_metadata <- merge(exome_bams_metadata, TSS, by='TSS_Code', all.x=TRUE)
  rm(TSS)

#7. Retrieve read depths of RNA-seq bams; read depths calculated by the script ***calc_read_depths_in_batches.txt***
  #Get metadata for all bams 
  RNAseq_bams_metadata <- read.table(paste('/home/mayo/m187735/s212975.Wickland_Immunomics/TCGA_metadata/RNAseq_bams/',CANCER_NAME,'_RNAseq_bams_metadata.txt',sep=''),header=TRUE)
  RNAseq_bams_metadata$submitter_id <- sub("^([^-]*-[^-]*-[^-]*-[^-]*).*", "\\1",RNAseq_bams_metadata$cases)
  
  #only interested in primary solid tumor
  RNAseq_bams_metadata  <- subset(RNAseq_bams_metadata, (tissue.definition=='Blood Derived Normal'))


#8. Combine dataframes
  #combine exome bams metadata with tumor expression data;keep all reads for subsequent merging with RNA-seq metadata
  combined_data_exome <- merge(exome_bams_metadata,read_depths,by='file_name')
  rm(exome_bams_metadata)
  #remove unnecessary columns
  combined_data_exome <- combined_data_exome[-c(1,2,3,4)]
  combined_data_exome <- combined_data_exome[,c(2,5,6,4,3,1)]
  #add suffix to columns
  names(combined_data_exome) <- c(paste(names(combined_data_exome)[1]),paste(names(combined_data_exome[,c(-1)]),'_exome',sep=''))
  
  #if any patient ID is duplicated, keep row/sample with highest read depth
  combined_data_exome_tmp <- aggregate(combined_data_exome[,2] ~ combined_data_exome[,1], combined_data_exome, max)
  names(combined_data_exome_tmp) <- c("submitter_id",names(combined_data_exome)[2])  
  combined_data_exome <- merge(combined_data_exome_tmp, combined_data_exome, by=c("submitter_id",names(combined_data_exome)[2]))
  rm (combined_data_exome_tmp)
  
  #combine RNAseq bams metadata with tumor expression data;keep all reads for subsequent merging with RNA-seq metadata
  combined_data_RNAseq <- merge(RNAseq_bams_metadata,read_depths,by='file_name')
  #remove unnecessary columns
  combined_data_RNAseq <- combined_data_RNAseq[-c(1,2,3,4)]
  #add suffix to columns
  names(combined_data_RNAseq) <- c(paste(names(combined_data_RNAseq)[1]),paste(names(combined_data_RNAseq[,c(-1)]),'_RNAseq',sep=''))
  rm(RNAseq_bams_metadata)
  rm(read_depths)
  
  #if any patient ID is duplicated, keep row/sample with highest read depth
  combined_data_RNAseq_tmp <- aggregate(combined_data_RNAseq[,2] ~ combined_data_RNAseq[,1], combined_data_RNAseq, max)
  names(combined_data_RNAseq_tmp) <- c("submitter_id",names(combined_data_RNAseq)[2])  
  combined_data_RNAseq <- merge(combined_data_RNAseq_tmp, combined_data_RNAseq, by=c("submitter_id",names(combined_data_RNAseq)[2]))
  rm (combined_data_RNAseq_tmp)
  
  #combine both bams+metadata dataframes
  combined_data <- merge(combined_data_exome,combined_data_RNAseq, by='submitter_id',all=TRUE)
  rm(combined_data_exome)
  rm(combined_data_RNAseq)
  

  #shorten case ID and remove any IDs that have duplicates
  combined_data$submitter_id_long <- combined_data$submitter_id
  combined_data$submitter_id <- sub("^([^-]*-[^-]*-[^-]*).*", "\\1",combined_data$submitter_id)
  
  combined_data <- combined_data[,c(1,3,4,5,6,7,8,2,9)]

  #merge with clinical data 
  combined_data <- merge(all_clinical_data,combined_data, by='submitter_id')
  rm(all_clinical_data)
  
  #merge with immune data
  combined_data <- merge(combined_data,immune_data, by='submitter_id',all.x=TRUE)
  rm(immune_data)


#9. Various modifications
  
  #clean up race, and keep only white / black / Asian
  combined_data$race[combined_data$race=="white"] <- "White"
  combined_data$race[combined_data$race=="black or african american"] <- "Black"
  combined_data$race[combined_data$race=="asian"] <- "Asian"
  combined_data$race[combined_data$race=="american indian or alaska native"] <- "Native"
  combined_data$race[combined_data$race=="not reported"] <- "Unknown"
  combined_data <- subset(combined_data, (race=="White") | (race=="Black") | (race=="Asian") | (race=="Unknown"))
  
  #consolidate tumor stage
  combined_data$tumor_stage <- as.character(combined_data$tumor_stage)
  combined_data$tumor_stage[combined_data$tumor_stage=="stage i"] <- "Stage I"
  combined_data$tumor_stage[combined_data$tumor_stage=="stage ia"] <- "Stage I"
  combined_data$tumor_stage[combined_data$tumor_stage=="stage ib"] <- "Stage I"
  combined_data$tumor_stage[combined_data$tumor_stage=="stage ii"] <- "Stage II"
  combined_data$tumor_stage[combined_data$tumor_stage=="stage iia"] <- "Stage II"
  combined_data$tumor_stage[combined_data$tumor_stage=="stage iib"] <- "Stage II"
  combined_data$tumor_stage[combined_data$tumor_stage=="stage iii"] <- "Stage III"
  combined_data$tumor_stage[combined_data$tumor_stage=="stage iiia"] <- "Stage III"
  combined_data$tumor_stage[combined_data$tumor_stage=="stage iiib"] <- "Stage III"
  combined_data$tumor_stage[combined_data$tumor_stage=="stage iiic"] <- "Stage III"
  combined_data$tumor_stage[combined_data$tumor_stage=="stage iv"] <- "Stage IV"
  combined_data$tumor_stage[combined_data$tumor_stage=="not reported"] <- NA
  combined_data$tumor_stage[combined_data$tumor_stage=="stage x"] <- NA
  
  
  #create log columns for TIL regional fraction, MANTIS score
  #combined_data$TIL_Regional_Fraction_log2 <- log(x=combined_data$TIL_Regional_Fraction,base=2)
  #combined_data$MANTIS_score_log2 <- log(x=combined_data$MANTIS_score, base=2)
  #combined_data$Mutation_Burden <- 2^combined_data$Mutation_Burden_log2
  
  #do below on mforge
  #FPKMs <- read.table('/research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/processing/TCGA/BRCA_expression/FPKM/all_genes_expression.txt',fill=TRUE)
  
  #plot(density(log(FPKMs[1:nrow(FPKMs),2])),ylim=c(0,0.13),main="Density plot of log2 FPKM from 1212 TCGA 'BRCA' samples")
  
  #for (i in 3:ncol(FPKMs)) {
  #  lines(density(log(FPKMs[1:nrow(FPKMs),i])))
  #}
  
  
  #add column for cancer
  combined_data$Cancer <- CANCER_NAME

#save dataset
write.table(x=combined_data,file=paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/overall_depths/',CANCER_NAME,'_overall_depths_etc_NORMAL.txt'),sep='\t',row.names=FALSE,quote=FALSE)
