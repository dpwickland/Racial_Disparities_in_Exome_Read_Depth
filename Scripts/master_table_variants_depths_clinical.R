#usage: Rscript --vanilla master_table_variants_depths_clinical.R <cancer name> <mapped reads or all reads>s
#in output, each line corresponds to ONE individual

options(warn=-1)

args<-commandArgs(TRUE)
.libPaths()
# test whether both arguments present
if (length(args)!=3) {
  stop("Usage example: Rscript --vanilla master_table_variants_depths_clinical.R BRCA mapped_reads RNAseq", call.=FALSE)
} else if (length(args)==3) {
  # default output file
  CANCER_NAME <- toString(args[1])
  MAPPED_OR_ALL <- toString(args[2])
  EXOME_OR_RNASEQ <- toString(args[3])
}

#load libraries
.libPaths(c( .libPaths(), "/home/mayo/m187735/R", "/usr/local/biotools/rpackages/R-3.5.2-latest", "/usr/local/biotools/rpackages/R-3.6.2-latest", "/usr/local/biotools/r/R-3.5.2/lib64/R/library"))
library(plyr)
library(dplyr, warn.conflicts=FALSE)
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

#3. Load tumor purity data (from Aran et al., 2015)    
tumor_purity_data <- data.frame(read.table('/home/mayo/m187735/s212975.Wickland_Immunomics/TCGA_metadata/tumor_purity_data_Aram.txt', header=TRUE))
#rename case ID column, tumor purity column
colnames(tumor_purity_data)[colnames(tumor_purity_data) == 'Sample_ID'] <- 'submitter_id'
colnames(tumor_purity_data)[colnames(tumor_purity_data) == 'CPE'] <- 'Tumor_Purity'
#subset cancer type
tumor_purity_data <- subset(tumor_purity_data, (Cancer_type == CANCER_NAME))
#remove superfluous column
tumor_purity_data <- tumor_purity_data[-c(2)]


#4. Load MSI (microsatellite) data from Bonneville
#MSI_Bonneville <- data.frame(read.table('/home/mayo/m187735/s212975.Wickland_Immunomics/TCGA_metadata#/Bonneville_MSI.txt', header=TRUE,sep='\t'))
#MSI_Bonneville <- subset(MSI_Bonneville, (Cancer.Type==paste("TCGA-",CANCER_NAME,sep='')))
#MSI_Bonneville <- MSI_Bonneville[-c(2)]
#names(MSI_Bonneville) <- c("submitter_id","MANTIS_score")


#5. Retrieve read depths of exome bams; read depths calculated by the script ***calc_exome_read_depths_in_batches.sh***
read_depths <- read.table(paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/overall_read_depths/',MAPPED_OR_ALL,'/',CANCER_NAME,'/all_samples_read_depths.txt'),header=FALSE)
colnames(read_depths) <- c("file_name","read_depth")
read_depths$bam_full_path <- read_depths$file_name 
read_depths$file_name <- sub(".*\\/", "\\1",read_depths$bam_full_path) #remove path from filename

#Get metadata for all bams 
bams_metadata <- read.table(paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/TCGA_metadata/',EXOME_OR_RNASEQ,'_bams/',CANCER_NAME,'_bams_metadata.txt'),header=TRUE)
bams_metadata$submitter_id <- sub("^([^-]*-[^-]*-[^-]*-[^-]*).*", "\\1",bams_metadata$cases)

#only interested in primary solid tumor or blood derived normal
bams_metadata  <- subset(bams_metadata, (tissue.definition=='Primary solid Tumor' | tissue.definition=='Blood Derived Normal'))

#extract tissue source site for exome
#extract center from barcode; must remove everything before last hyphen
bams_metadata$center <- sub(".*\\-", "\\1",bams_metadata$cases)

#extract TSS code **from barcode**; remove everything after 2nd hyphen
bams_metadata$TSS_Code <- gsub("TCGA-","",sub("^([^-]*-[^-]*).*", "\\1",bams_metadata$cases))

#match TSS code with site ID (from https://gdc.cancer.gov/files/public/file/tcga_code_tables.zip) -- replace spaces, apostrophes with underscores first
TSS <- read.table('/home/mayo/m187735/s212975.Wickland_Immunomics/TCGA_metadata/tissue_source_site.txt', header=TRUE)
bams_metadata <- merge(bams_metadata, TSS, by='TSS_Code', all.x=TRUE)
rm(TSS)

#match center code with center ID (same source as above)
SS <- read.table('/home/mayo/m187735/s212975.Wickland_Immunomics/TCGA_metadata/tissue_seq_center.txt',header=TRUE)
bams_metadata$center <- as.numeric(bams_metadata$center)
bams_metadata <- merge(bams_metadata,SS, by='center',all.x=TRUE)

#7. Retrieve read depths of RNA-seq bams; read depths calculated by the script ***calc_read_depths_in_batches.txt***
#Get metadata for all bams 
#  RNAseq_bams_metadata <- read.table(paste('/home/mayo/m187735/s212975.Wickland_Immunomics/TCGA_metadata/RNAseq_bams/',CANCER_NAME,'_RNAseq_bams_metadata.txt',sep=''),header=TRUE)
#  RNAseq_bams_metadata$submitter_id <- sub("^([^-]*-[^-]*-[^-]*-[^-]*).*", "\\1",RNAseq_bams_metadata$cases)

#only interested in primary solid tumor
#  RNAseq_bams_metadata  <- subset(RNAseq_bams_metadata, (tissue.definition=='Primary solid Tumor' | tissue.definition=='Blood_Derived_Normal'))


#8. Combine dataframes
#combine exome bams metadata with tumor expression data;keep all reads for subsequent merging with RNA-seq metadata
combined_data <- merge(bams_metadata,read_depths,by='file_name')
rm(bams_metadata)
#remove unnecessary columns
combined_data <- combined_data[-c(1,2,3,4,5)]
combined_data <- combined_data[,c(2,5,6,4,3,1)]
#add suffix to columns
#names(combined_data) <- c(paste(names(combined_data)[1]),paste(names(combined_data[,c(-1)]),'_',EXOME_OR_RNASEQ,sep=''))

#if any patient ID is duplicated, keep row/sample with highest read depth
combined_data_tmp <- aggregate(combined_data[,2] ~ combined_data[,1], combined_data, max)
names(combined_data_tmp) <- c("submitter_id",names(combined_data)[2])  
combined_data <- merge(combined_data_tmp, combined_data, by=c("submitter_id",names(combined_data)[2]))
rm (combined_data_tmp)

#combine RNAseq bams metadata with tumor expression data;keep all reads for subsequent merging with RNA-seq metadata
#  combined_data_RNAseq <- merge(RNAseq_bams_metadata,read_depths,by='file_name')
#remove unnecessary columns
#  combined_data_RNAseq <- combined_data_RNAseq[-c(1,2,3,4)]
#add suffix to columns
#  names(combined_data_RNAseq) <- c(paste(names(combined_data_RNAseq)[1]),paste(names(combined_data_RNAseq[,c(-1)]),'_RNAseq',sep=''))
#  rm(RNAseq_bams_metadata)
#  rm(read_depths)

#if any patient ID is duplicated, keep row/sample with highest read depth
#  combined_data_RNAseq_tmp <- aggregate(combined_data_RNAseq[,2] ~ combined_data_RNAseq[,1], combined_data_RNAseq, max)
#  names(combined_data_RNAseq_tmp) <- c("submitter_id",names(combined_data_RNAseq)[2])  
#  combined_data_RNAseq <- merge(combined_data_RNAseq_tmp, combined_data_RNAseq, by=c("submitter_id",names(combined_data_RNAseq)[2]))
#  rm (combined_data_RNAseq_tmp)

#combine both bams+metadata dataframes
#combined_data <- merge(combined_data,combined_data_RNAseq, by='submitter_id',all=TRUE)

# rm(combined_data_RNAseq)

#combined with tumor purity data
combined_data <- merge(combined_data, tumor_purity_data, by='submitter_id',all.x=TRUE)
rm(tumor_purity_data)

#shorten case ID and remove any IDs that have duplicates
combined_data$submitter_id_long <- combined_data$submitter_id
combined_data$submitter_id <- sub("^([^-]*-[^-]*-[^-]*).*", "\\1",combined_data$submitter_id)


#if any patient ID is duplicated, keep row/sample with highest tumor purity value (most pairs have just one member with a tumor purity value anyway)
combined_data_tmp <- aggregate(combined_data[,7] ~ combined_data[,8], combined_data, max)
names(combined_data_tmp) <- c("submitter_id_long",names(combined_data)[7])  
combined_data <- merge(combined_data_tmp, combined_data, by=c("submitter_id_long",names(combined_data)[7]),all.y=TRUE)
combined_data <- combined_data[,c(1,3,4,5,6,7,2,8)]
rm(combined_data_tmp)

#also, barcode must end in 01A or 10A
combined_data <- combined_data %>% filter(str_detect(combined_data$submitter_id_long, "-01A|-10A"))

#merge with MSI data
#combined_data <- merge(combined_data, MSI_Bonneville, by='submitter_id')
#rm(MSI_Bonneville)

#merge with clinical data 
combined_data <- merge(combined_data,all_clinical_data, by='submitter_id',all.x=TRUE)
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

#tissue type
combined_data$tissue.definition<- sub("Blood Derived Normal","Blood_Normal",combined_data$tissue.definition)
combined_data$tissue.definition <- sub("Primary solid Tumor","Primary_Tumor",combined_data$tissue.definition)


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

#add capture kit  
# data from here: https://api.gdc.cancer.gov/files?filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22!=%22%2C%22content%22%3A%7B%22field%22%3A%22cases.submitter_id%22%2C%22value%22%3A%5B%22TCGA-4V-MIUM%22%5D%7D%7D%2C%7B%22op%22%3A%22%3D%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_type%22%2C%22value%22%3A%22Aligned%20Reads%22%7D%7D%5D%7D&format=tsv&fields=file_id,file_name,cases.submitter_id,cases.case_id,data_category,data_type,cases.samples.tumor_descriptor,cases.samples.tissue_type,cases.samples.sample_type,cases.samples.submitter_id,cases.samples.sample_id,analysis.workflow_type,cases.project.project_id,cases.samples.portions.analytes.aliquots.aliquot_id,cases.samples.portions.analytes.aliquots.submitter_id,analysis.metadata.read_groups.target_capture_kit_name,analysis.metadata.read_groups.target_capture_kit_target_region&size=100000
#these date are consistent with the BUckley 2017 paper, which howevr lists only the normal samples

capture_kit <- read.csv('/home/mayo/m187735/s212975.Wickland_Immunomics/TCGA_metadata/capture_kit_TCGA_bamz.csv',header=TRUE)

combined_data <- merge(combined_data, capture_kit[,c(which(colnames(capture_kit)=='SUBMITTER_ID_MED' | colnames(capture_kit)=='CAPTURE_KIT'))], by.x='submitter_id_long',by.y='SUBMITTER_ID_MED')

#combined_data <- unique(combined_data)
#combined_data$capture_kit_used <- 'Yes'

#now turn that capture kit column into multiple
#combined_data <- combined_data %>% spread(CAPTURE_KIT,capture_kit_used)

#save dataset
write.table(x=combined_data,file=paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/overall_read_depths_and_sample_information/',CANCER_NAME,'_',MAPPED_OR_ALL,'_',EXOME_OR_RNASEQ,'_master_table.txt'),sep='\t',row.names=FALSE,quote=FALSE)