#usage: Rscript --vanilla TCGA_gene_depths_pairs_from_VCF.R BRCA TP53

args<-commandArgs(TRUE)
.libPaths()
# test whether both arguments present
if (length(args)!=1) {
  stop("Usage example: Rscript --vanilla TCGA_gene_depths_pairs_from_VCF.R TP53", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  GENE <- toString(args[1])
}

#load libraries
.libPaths(c( .libPaths(), "/home/mayo/m187735/R", "/usr/local/biotools/rpackages/R-3.5.2-latest",  "/usr/local/biotools/r/R-3.5.2/lib64/R/library"))
library(plyr)
library(dplyr)
library(readr)
library(stringr)
library(gtools)
setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)

############################################
#INDIVIDUAL GENE ANALYSIS: EACH ROW IS A DETECTED SOMATIC MUTATION FROM ONE PARTICULAR PATIENT#
############################################
#0. Retrieve gene coordinates from gtf file (downloaded from wget https://api.gdc.cancer.gov/data/25aa497c-e615-4cb7-8751-71f744f9691f; then paste <(zcat 25aa497c-e615-4cb7-8751-71f744f9691f | awk '$3 == "gene"' | cut -f1,4,5,7)  <(zcat 25aa497c-e615-4cb7-8751-71f744f9691f | awk '$3 == "gene"' | cut -d ' ' -f2,8)  > TCGA_GRCH38_v22.gtf
gtf <- read.table('/home/mayo/m187735/s212975.Wickland_Immunomics/TCGA_GRCH38_v22.gtf',header=FALSE)
#remove blank cols
gtf <- gtf[-c(6,8)]
#name cols
colnames(gtf) <- c("chr","start","end","strand","ensemble_ID","gene_name")
#remove trailing part of ensemble ID
gtf$ensemble_ID <- sub("^([^.]*).*", "\\1",gtf$ensemble_ID)

#subset gtf for gene of interest
gene_subset <- subset(gtf, (gene_name==GENE))
rm(gtf)
GENE_AND_ENSEMBL_ID <- c(gene_subset[,5])
names(GENE_AND_ENSEMBL_ID) <- GENE
DRIVER_GENES <- c(c(as.character(gene_subset[,1]),gene_subset[,2],gene_subset[,3]))
names(DRIVER_GENES) <- c(paste(GENE,"_1",sep=''),paste(GENE,"_2",sep=''),paste(GENE,"_3",sep=''))

CANCER_LIST <- c('BRCA','LUAD','UCEC','COAD','PRAD','KIRC')

for (CANCER_NAME in CANCER_LIST){
    
  #1. Extract mutations in particular genes from somatic VCFs (which contain the ensemble ID)
    system(paste("mkdir -p /home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/mutation_lists/",CANCER_NAME,"/",sub("^([^_]*).*", "\\1",names(DRIVER_GENES)[1]),sep=''))  
    system(paste("for file in /home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/VCFs/all_mutec_VCFs/",CANCER_NAME,"/*/*.vcf.gz; do echo $file; zcat $file | awk '$1==\"",DRIVER_GENES[1],"\" && $2 >=",DRIVER_GENES[2]," && $2 <=",DRIVER_GENES[3],"' | cut -f1,2,4,5,11 > /home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/mutation_lists/",CANCER_NAME,"/",sub("^([^_]*).*", "\\1",names(DRIVER_GENES)[1]),"/`basename ${file} .vcf.gz`_",sub("^([^_]*).*", "\\1",names(DRIVER_GENES)[1]),"_mutations.txt; done;",sep='')) 
  
  #2. Aggregate mutations found above into single file, adding full path to the files
    system(paste("for file in /home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/mutation_lists/",CANCER_NAME,"/",sub("^([^_]*).*", "\\1",names(DRIVER_GENES)[1]),"/*vep_",sub("^([^_]*).*", "\\1",names(DRIVER_GENES)[1]),"_mutations.txt; do awk '{print (FILENAME (NF?\"\t\":\"\")) $0}' $file | sed 's/.vep_",sub('^([^_]*).*', '\\1',names(DRIVER_GENES)[1]),"_mutations.txt/.vep.vcf.gz/g';done > /home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/mutation_lists/",CANCER_NAME,"/all_",sub('^([^_]*).*', '\\1',names(DRIVER_GENES)[1]),"_mutations.txt",sep=''))
    
    cdriver_tumor_mutations <- read.table(paste("/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/mutation_lists/",CANCER_NAME,"/all_",sub('^([^_]*).*', '\\1',names(DRIVER_GENES)[1]),"_mutations.txt",sep=''),header=FALSE)
    
    #only keep SNVs -- no INDELs
    cdriver_tumor_mutations <- subset(cdriver_tumor_mutations, (str_length(cdriver_tumor_mutations$V4)==1) & (str_length(cdriver_tumor_mutations$V5)==1))
    
    #create copy of this dataframe for use in an analysis of all pairwise comparisons between individual bams and the list of somatic mutations detected in the VCFs 
    #system(paste0("mkdir /home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/gene_depths_pairs_from_VCF/",GENE))
    #write.table(x=cdriver_tumor_mutations,file=paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/gene_depths_pairs_from_VCF/',GENE,'/',CANCER_NAME,'_',GENE,'_all_detected_mutations.txt'),sep='\t')
    
    
    #remove path to VCFs mutation counts files
    cdriver_tumor_mutations$V1 <- gsub(paste("/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/mutation_lists/",CANCER_NAME,"/",sub('^([^_]*).*', '\\1',names(DRIVER_GENES)[1]),"/",sep=''),"",cdriver_tumor_mutations$V1)
    
    #isolate genotype from VCF
    cdriver_tumor_mutations$V6 <- sub("^([^:]*:[^:]*).*", "\\1",cdriver_tumor_mutations$V6)
    
    #isolate allelic depths from VCF
    cdriver_tumor_mutations$GTs <- sub("^([^:]*).*", "\\1",cdriver_tumor_mutations$V6) #delete everything after second colon
    cdriver_tumor_mutations$ADs <- sub(".*\\:", "\\1",cdriver_tumor_mutations$V6) #delete everything up to and including first colon
    cdriver_tumor_mutations$V8 <- sub("^([^,]*).*", "\\1",cdriver_tumor_mutations$ADs) #delete everything after comma
    cdriver_tumor_mutations$V9 <- sub(".*\\,", "\\1",cdriver_tumor_mutations$ADs) #delete everything before comma
    
    #rename, remove cols
    cdriver_tumor_mutations <- cdriver_tumor_mutations[-c(6,8)] #remove AD col
    names(cdriver_tumor_mutations) <- c("file_name","chromosome","mutation_position","ref_allele_VCF","alt_allele_VCF","genotype_VCF","AD_ref_VCF","AD_alt_VCF")
    
    #add additional columns: total AD, alt allele concentration (from AD)
    cdriver_tumor_mutations$AD_total_VCF <- as.numeric(cdriver_tumor_mutations$AD_ref_VCF) + as.numeric(cdriver_tumor_mutations$AD_alt_VCF)
    cdriver_tumor_mutations$alt_allele_concent_VCF <- as.numeric(cdriver_tumor_mutations$AD_alt_VCF) / as.numeric(cdriver_tumor_mutations$AD_total_VCF)
  
  
  #3. Match with TCGA metadata
    all_VCF_metadata <- read.table(paste('/home/mayo/m187735/s212975.Wickland_Immunomics/TCGA_metadata/somatic_VCF/',CANCER_NAME,'_VCF_metadata.txt',sep=''),header=TRUE)
    all_VCF_metadata <- (subset(all_VCF_metadata,(data_type=="Annotated Somatic Mutation" & analysis_workflow_type=="MuTect2 Annotation")))
    all_VCF_metadata <- all_VCF_metadata[,c("file_name","submitter_id")]
    
    #get tissue definition from barcodes b/c one from metadata is wrong
    tissue1 <- sub("_.*", "", all_VCF_metadata$submitter_id) #get first bacode
    tissue1 <- sub("^([^-]*-[^-]*-[^-]*-[^-]*).*", "\\1",tissue1) #delete everything after and including 4th hyphen
    tissue1 <- sub(".*\\-", "\\1",tissue1) #delete everything before last hyphen
    tissue1 <- sub('.{1}$', '', tissue1) #delete last character
    tissue1 <- gsub("01","Primary solid Tumor",tissue1)
    tissue1 <- gsub("06","Metastatic",tissue1)
    
    tissue2 <- gsub("_mutect","",gsub("^[^_]+_|_[^_]+$", "", all_VCF_metadata$submitter_id)) #get second barcode
    tissue2 <- sub("^([^-]*-[^-]*-[^-]*-[^-]*).*", "\\1",tissue2) #delete everything after and including 4th hyphen
    tissue2 <- sub(".*\\-", "\\1",tissue2) #delete everything before last hyphen
    tissue2 <- sub('.{1}$', '', tissue2) #delete last character
    tissue2 <- gsub("10","Blood Derived Normal",tissue2)
    tissue2 <- gsub("11","Solid Tissue Normal",tissue2)
    
    #assign tissue to VCF metadata dataframe
    all_VCF_metadata$tissue1 <- tissue1
    all_VCF_metadata$tissue2 <- tissue2
    rm(tissue1)
    rm(tissue2)
    
    #rename cols
    names(all_VCF_metadata) <- c(paste(names(all_VCF_metadata)[1:2]),paste("somatic_VCF",names(all_VCF_metadata[,c(-1,-2)]),sep='_'))
    
    #shorten submitter id - will contain only the primary tumor barcode or metastatic barcode
    all_VCF_metadata$submitter_id_long <- sub("^([^-]*-[^-]*-[^-]*-[^-]*).*", "\\1",all_VCF_metadata$submitter_id)
    all_VCF_metadata$submitter_id <- sub("^([^-]*-[^-]*-[^-]*).*", "\\1",all_VCF_metadata$submitter_id)
    all_VCF_metadata <- subset(all_VCF_metadata, (somatic_VCF_tissue1=='Primary solid Tumor'))
    
    #merge mutations dataframe with VCF metadata dataframe; note that not all patients in cdriver_tumor_mutations have VCF metadata
    #DOUBLE CHECK THIS -- WHY DO THESE FILES HAVE NO VCF METADATA?
    cdriver <- merge(cdriver_tumor_mutations, all_VCF_metadata, by='file_name')
    
    #delete extra dataframes
    rm(cdriver_tumor_mutations)
    rm(all_VCF_metadata)
    
    #some patients have multiple VCFs; in such instances, keep only the VCF with greatest read depth
    cdriver_tmp1 <- aggregate(cdriver$AD_total_VCF ~ cdriver$file_name, cdriver, max)
    names(cdriver_tmp1) <- c('file_name','AD_total_VCF')
    cdriver <- merge(cdriver_tmp1, cdriver, by.x=c('file_name','AD_total_VCF'))
    
    #some will have same read depths in mult. files -- in those cases, just keep the first one listed
    cdriver <- cdriver[!duplicated(cdriver[,c('submitter_id','mutation_position')]),]
  
    #remove file_name col
    cdriver <- cdriver[-c(1)]
    
    #place submitter IDs first
    cdriver <- cdriver[,c(10,ncol(cdriver),2:8,1,9,12:ncol(cdriver)-1)]
    
  #4. Add exome-bams-metadata and exome bam paths (from overall read depths data generated by previous script)
  
    #Get metadata for all bams; duplicates fine b/c some patients will have multiple mutations 
    exome_bams_metadata <- read.table(paste('/home/mayo/m187735/s212975.Wickland_Immunomics/TCGA_metadata/exome_bams/',CANCER_NAME,'_bams_metadata.txt',sep=''),header=TRUE)
    exome_bams_metadata$submitter_id <- sub("^([^-]*-[^-]*-[^-]*-[^-]*).*", "\\1",exome_bams_metadata$cases)
    
    #only interested in primary solid tumor
    exome_bams_metadata  <- subset(exome_bams_metadata, (tissue.definition=='Primary solid Tumor'))
    
    colnames(exome_bams_metadata)[colnames(exome_bams_metadata) == 'submitter_id'] <- 'submitter_id_long'
    cdriver <- merge(cdriver,exome_bams_metadata,by='submitter_id_long',all.x=TRUE)
    rm(exome_bams_metadata)
    
    #subset only the primary solid tumor bams --i.e.the indivs from VCFs for whom we have bams
    cdriver <- subset(cdriver, (tissue.definition=='Primary solid Tumor'))
    
    
    #Get overall read depths data
    combined_data <- read.table(paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/overall_depths/',CANCER_NAME,'_overall_depths_etc.txt'),header=TRUE)
    
    #perhaps better to combine based on exome bam filename, rather than on submitter_id????
    combined_data$exome_filename_without_path <- sub(".*\\/", "\\1",combined_data$bam_full_path_exome) #delete path
    
    #combine combined_data (containing overall read depths) and cdriver data; some patients with VCFs may not have bams, and vice versa; combined based on exome file name, rather than submitter_id (less ambiguity)
    cdriver <- merge(cdriver, combined_data, by.x='file_name',by.y='exome_filename_without_path')  
  
    
    #if any patient ID has multiple exome files, keep only the read depth file with greatest depth
    cdriver_tmp2 <- aggregate(cdriver$read_depth_exome ~ cdriver$submitter_id.x, cdriver, max)
    names(cdriver_tmp2) <- c("submitter_id.x","read_depth_exome")
    cdriver <- merge(cdriver_tmp2, cdriver, by=c('submitter_id.x','read_depth_exome'))
    names(cdriver)[1] <- 'submitter_id'
  
    #add col for number of mutations per sample
    cdriver <- add_count(cdriver, submitter_id)
    names(cdriver)[ncol(cdriver)] <- 'mutation_count_for_patient'
    
    
  #5. Read depth calculations: Count number of mapped reads of each cdriver mutated position in each individual
    variant_read_depth_bam_VIEW_col <- ncol(cdriver)+1
    variant_read_depth_bam_MPILEUP_col <- ncol(cdriver)+2
    ref_base_count_bam_col <- ncol(cdriver)+3
    alt_base_count_bam_col <- ncol(cdriver)+4
    #avg_read_mapping_quality_ref_bam_col <- ncol(cdriver)+5
    #avg_read_mapping_quality_alt_bam_col <- ncol(cdriver)+6
    avg_base_quality_ref_bam_col <- ncol(cdriver)+5
    avg_base_quality_alt_bam_col <- ncol(cdriver)+6
  
    for (i in 1:nrow(cdriver)){ 
      print(paste0(CANCER_NAME,": ",i,"/",nrow(cdriver)))
      #total_read_depth <- numeric(0) #initialize object to store depths for particular case ID / row
      #row_position_read_depth <- numeric(0) #initialize object to store depths for particular case ID / row
      #output NA if index file does not exist
      variant_read_depth_bam_VIEW <- system(paste("if [ -f ",sub(".bam",".bai",cdriver$bam_full_path_exome[i])," ]; then /research/bsi/tools/biotools/samtools/1.9/miniconda/bin/samtools view -c -F4 -F1024 -q 20",cdriver$bam_full_path_exome[i],paste(cdriver$chromosome[i],":",cdriver$mutation_position[i],"-",cdriver$mutation_position[i],"; else echo NA; fi",sep='')),intern=TRUE,ignore.stderr=TRUE)
      cdriver[i,variant_read_depth_bam_VIEW_col] <- variant_read_depth_bam_VIEW
      
      pileup <- system(paste("if [ -f ",sub(".bam",".bai",cdriver$bam_full_path_exome[i])," ]; then /research/bsi/tools/biotools/samtools/1.9/miniconda/bin/samtools mpileup ",cdriver$bam_full_path_exome[i],"--min-MQ 20 --min-BQ 0 --count-orphans --output-MQ --region ",paste(cdriver$chromosome[i],":",cdriver$mutation_position[i],"-",cdriver$mutation_position[i],"; else echo NA; fi",sep='')),intern=TRUE,ignore.stderr=TRUE)
      
      if (length(pileup)!=0 && pileup!= 'NA'){
        variant_read_depth_bam_MPILEUP <- strsplit(pileup,'\t')[[1]][4] 
        bases <- strsplit(pileup,'\t')[[1]][5]
        bases <- gsub("[^tcgaTCGA+]+","",bases) #remove $ character from list of bases; this simply; don't need to do this for the qualities b/c no. qual values is already equal to no. bases
  
        ref_base_count_bam <- str_count(toupper(bases),toupper(as.character(cdriver[i,]$ref_allele_VCF))) #count number of reads with ref base in the pileup (convert all to uppercase)
        alt_base_count_bam <- str_count(toupper(bases),toupper(as.character(cdriver[i,]$alt_allele_VCF))) #count number of reads with alt base in the pileup (convert all to uppercase)
       
        cdriver[i,variant_read_depth_bam_MPILEUP_col] <- variant_read_depth_bam_MPILEUP
        cdriver[i,ref_base_count_bam_col] <- ref_base_count_bam
        cdriver[i,alt_base_count_bam_col] <- alt_base_count_bam
  
        ref_indices <- unlist(str_locate_all(toupper(bases),toupper(as.character(cdriver[i,]$ref_allele_VCF))))[1:ref_base_count_bam] #identify index positions of REF allele in the pileup; then get their corresponding quality scores; unlisting must be chopped in half b/c entire vector repeated once
        alt_indices <- unlist(str_locate_all(toupper(bases),toupper(as.character(cdriver[i,]$alt_allele_VCF))))[1:alt_base_count_bam] #identify index positions of REF allele in the pileup; then get their corresponding quality scores; unlisting must be chopped in half b/c entire vector repeated once
      
    #    read_mapping_quality <- asc(strsplit(pileup,'\t')[[1]][7])-33
    #    avg_read_mapping_quality_ref_bam <- round(mean(read_mapping_quality[ref_indices]),digits=2)
    #    avg_read_mapping_quality_alt_bam <- round(mean(read_mapping_quality[alt_indices]),digits=2) #problem: insertion sequences; see http://samtools.sourceforge.net/pileup.shtml; so for now, just don't do mapping quality
        
        base_quality <- asc(strsplit(pileup,'\t')[[1]][6])-33
        avg_base_quality_ref_bam <- round(mean(base_quality[ref_indices]),digits=2)
        avg_base_quality_alt_bam <- round(mean(base_quality[alt_indices]),digits=2)
        
  
        
     #  cdriver[i,avg_read_mapping_quality_ref_bam_col] <- avg_read_mapping_quality_ref_bam
     #  cdriver[i,avg_read_mapping_quality_alt_bam_col] <- avg_read_mapping_quality_alt_bam
        
        cdriver[i,avg_base_quality_ref_bam_col] <- avg_base_quality_ref_bam
        cdriver[i,avg_base_quality_alt_bam_col] <- avg_base_quality_alt_bam
        
      }
      
      else if (length(pileup)==0 || pileup=='NA'){
        cdriver[i,variant_read_depth_bam_MPILEUP_col] <- NA
        cdriver[i,ref_base_count_bam_col] <- NA
        cdriver[i,alt_base_count_bam_col] <- NA
        
        
        #cdriver[i,avg_read_mapping_quality_ref_bam_col] <- NA
        #cdriver[i,avg_read_mapping_quality_alt_bam_col] <- NA
        cdriver[i,avg_base_quality_ref_bam_col] <- NA
        cdriver[i,avg_base_quality_alt_bam_col] <- NA
        
      }
      
      #total_read_depth <- system(paste("/research/bsi/tools/biotools/samtools/1.9/miniconda/bin/samtools view -c -F4",cdriver$file_name_exome_bams[i]),intern=TRUE)
      #cdriver[i,total_read_depth_col] <- total_read_depth
      i <- i + 1 #increase row
    }
    
    colnames(cdriver)[variant_read_depth_bam_VIEW_col] <- 'variant_read_depth_bam_VIEW'
    colnames(cdriver)[variant_read_depth_bam_MPILEUP_col] <- 'variant_read_depth_bam_MPILEUP'
    colnames(cdriver)[ref_base_count_bam_col] <- 'ref_base_count_bam'
    colnames(cdriver)[alt_base_count_bam_col] <- 'alt_base_count_bam'
    #colnames(cdriver)[avg_read_mapping_quality_ref_bam_col] <- 'avg_read_mapping_quality_ref_bam'
    #colnames(cdriver)[avg_read_mapping_quality_alt_bam_col] <- 'avg_read_mapping_quality_alt_bam'
    colnames(cdriver)[avg_base_quality_ref_bam_col] <- 'avg_ref_base_quality_bam'
    colnames(cdriver)[avg_base_quality_alt_bam_col] <- 'avg_alt_base_quality_bam'
    
    #add additional column: alt allele concentration from bams
    cdriver$alt_allele_concent_bam <- as.numeric(cdriver$alt_base_count_bam) / as.numeric(cdriver$variant_read_depth_bam_MPILEUP)  
    
    #remove extra columns (submitter_id_long.x,file_name,file_id,cases,cases,submitter_id.y,submitter_id_long.y); then rearrange
    cdriver <- cdriver[-c(4,3,16,17,19,28)]
    cdriver <- cdriver[,c(1,3:10,2,11:ncol(cdriver))]
  
    #save dataset
    system(paste0('mkdir /home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/gene_depths_pairs_from_VCF/',GENE,'/'))
    system(paste0('rm /home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/gene_depths_pairs_from_VCF/',GENE,'/',CANCER_NAME,'_',GENE,'_depths_etc.txt'))
    write.table(x=cdriver,file=paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/gene_depths_pairs_from_VCF/',GENE,'/',CANCER_NAME,'_',GENE,'_depths_etc.txt'),sep='\t',row.names=FALSE)
    
    
  #6. Also output FPKM for each patient, both tumor and normal
    
    for (GENE in GENE_AND_ENSEMBL_ID) {
      
      #make new directory for gene (delete directory if already exists)
      system(paste("rm -r /home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/FPKM/",CANCER_NAME,"/",names(GENE_AND_ENSEMBL_ID)[GENE_AND_ENSEMBL_ID==GENE],"_",GENE,sep=''))
      system(paste("mkdir -p /home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/FPKM/",CANCER_NAME,"/",names(GENE_AND_ENSEMBL_ID)[GENE_AND_ENSEMBL_ID==GENE],"_",GENE,sep=''))
      
      #retrieve expression levels from each sample for gene of interest (from TCGA FPKM files)
      system(paste("for file in /home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/FPKM/",CANCER_NAME,"/input/*/*gz; do zcat $file | grep ",GENE," > /home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/FPKM/",CANCER_NAME,"/",names(GENE_AND_ENSEMBL_ID)[GENE_AND_ENSEMBL_ID==GENE],"_",GENE,"/`basename ${file} .FPKM.txt.gz`;done",sep=''))
      
      #for each gene, output file listing file name and expression for that gene
      system(paste("for file in /home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/FPKM/",CANCER_NAME,"/",names(GENE_AND_ENSEMBL_ID)[GENE_AND_ENSEMBL_ID==GENE],"_",GENE,"/*; do awk '{print basename (FILENAME (NF?\"\\t\":\"\")) $0}' ${file} | cut -f1,3 > /${file}_trimmed;done",sep=''))
      
      #aggregate info from all patients and all genes
      gene_files = list.files(path=paste("/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/FPKM/",CANCER_NAME,"/",names(GENE_AND_ENSEMBL_ID)[GENE_AND_ENSEMBL_ID==GENE],"_",GENE,"/",sep=''), pattern="*_trimmed", full.names=TRUE)
      all_patients_gene_FPKM = ldply(gene_files, read.table, sep = "\t", fill=TRUE, header = FALSE)
      
      #make col names, deleting all parts of file path
      colnames(all_patients_gene_FPKM) <- c("file_name",paste(names(GENE_AND_ENSEMBL_ID)[GENE_AND_ENSEMBL_ID==GENE],"_",GENE,sep=''))
      all_patients_gene_FPKM$file_name <- sub(".*\\/", "\\1",all_patients_gene_FPKM$file_name)
      
      if (GENE==GENE_AND_ENSEMBL_ID[1]) {
        selected_genes_expression_all_patients <- all_patients_gene_FPKM
      }
      else if (GENE!=GENE_AND_ENSEMBL_ID[1]) {
        selected_genes_expression_all_patients <- merge(all_patients_gene_FPKM,selected_genes_expression_all_patients,by='file_name')
      }
    }
    
    rm(all_patients_gene_FPKM)
    expression_data <- selected_genes_expression_all_patients
    rm(selected_genes_expression_all_patients)
    
    #get metadata for the FPKM files retrieved above
    TCGA_expression_metadata <- read.table(paste('/home/mayo/m187735/s212975.Wickland_Immunomics/TCGA_metadata/expression_FPKM/',CANCER_NAME,'_expression_metadata.txt',sep=''),header=TRUE)
    
    #extract key information from sample barcodes in metadata object
    TCGA_expression_metadata <- TCGA_expression_metadata[,c("file_name","cases","tissue.definition")]
    
    #extract center from barcode; must remove everything before last hyphen
    TCGA_expression_metadata$center <- sub(".*\\-", "\\1",TCGA_expression_metadata$cases)
    
    #extract TSS code **from barcode**; remove everything after 2nd hyphen
    TCGA_expression_metadata$TSS_Code <- gsub("TCGA-","",sub("^([^-]*-[^-]*).*", "\\1",TCGA_expression_metadata$cases))
    
    #match TSS code with site ID (from https://gdc.cancer.gov/files/public/file/tcga_code_tables.zip) -- replace spaces, apostrophes with underscores first
    TSS <- read.table('/home/mayo/m187735/s212975.Wickland_Immunomics/TCGA_metadata/tissue_source_site.txt', header=TRUE)
    TCGA_expression_metadata <- merge(TCGA_expression_metadata, TSS, by='TSS_Code', all.x=TRUE)
    
    #delete superflous cols in TCGA_expression_metadata
    TCGA_expression_metadata <- TCGA_expression_metadata[-c(1,5)]
    
    #shorten case ID; remove everything after 4th hyphen (to integrate with clinical)
    TCGA_expression_metadata$cases <- sub("^([^-]*-[^-]*-[^-]*-[^-]*).*", "\\1",TCGA_expression_metadata$cases) 
    
    #remove filename extension and rename columns
    TCGA_expression_metadata$file_name <- sub("^([^.]*).*", "\\1",TCGA_expression_metadata$file_name) 
    colnames(TCGA_expression_metadata)[colnames(TCGA_expression_metadata) == 'cases'] <- 'submitter_id'
    
    #combine above dataframes on file_name, remove filename, and delete uncombined
    all_expression <- merge(TCGA_expression_metadata,expression_data,  by='file_name')
    all_expression <- all_expression[-c(1)]
    rm(TCGA_expression_metadata)
    rm(expression_data)
    
    # create separate dataframes for tumor, normal, with special suffixes for each
    #tumor
    tumor_expression <- subset(all_expression, (tissue.definition == "Primary solid Tumor"))
    names(tumor_expression) <- c(paste(names(tumor_expression)[1]),paste(names(tumor_expression[,c(-1)]),'_FPKM_tumor',sep=''))
    
      #if any patient ID is duplicated, keep row/sample with highest FPKM value
      tumor_expression_tmp <- aggregate(tumor_expression[,4] ~ tumor_expression[,1], tumor_expression, max)
      names(tumor_expression_tmp) <- c("submitter_id",names(tumor_expression)[4])  
      tumor_expression <- merge(tumor_expression_tmp, tumor_expression, by=c("submitter_id",names(tumor_expression)[4]))
      tumor_expression <- tumor_expression[,c(1,3,4,2)]
      rm (tumor_expression_tmp)
    
    #normal
    normal_expression <- subset(all_expression, (tissue.definition == "Solid Tissue Normal"))
    names(normal_expression) <- c(paste(names(normal_expression)[1]),paste(names(normal_expression[,c(-1)]),'_FPKM_normal',sep=''))
    
      #if any patient ID is duplicated, keep row/sample with highest FPKM value
      normal_expression_tmp <- aggregate(normal_expression[,4] ~ normal_expression[,1], normal_expression, max)
      names(normal_expression_tmp) <- c("submitter_id",names(normal_expression)[4])  
      normal_expression <- merge(normal_expression_tmp, normal_expression, by=c("submitter_id",names(normal_expression)[4]))
      normal_expression <- normal_expression[,c(1,3,4,2)]
      rm(normal_expression_tmp)
    
    #shorten case ID for normal
    normal_expression$submitter_id <- sub("^([^-]*-[^-]*-[^-]*).*", "\\1",normal_expression$submitter_id)
    rm(all_expression) 
    
  
    #merge with tumor expression with overall depths dataset
    colnames(tumor_expression)[colnames(tumor_expression) == 'submitter_id'] <- 'submitter_id_long'
    combined_data <- merge(combined_data, tumor_expression, by='submitter_id_long',all.x=TRUE)
    rm(tumor_expression)
    
    #merge with normal expression; keep all from combined
    combined_data <- merge(combined_data, normal_expression,by='submitter_id',all.x=TRUE)
    rm(normal_expression)   
    
    #based on above, throw out values below 0 when log taken (but don't keep log-transformed values...yet)
  #  for (GENE in GENE_AND_ENSEMBL_ID) {
      #on tumor samples
    #  gene_col_tumor <- paste(names(GENE_AND_ENSEMBL_ID)[GENE_AND_ENSEMBL_ID==GENE],"_",GENE,"_FPKM_tumor",sep="")  
    #  combined_data[,c(gene_col_tumor)][log(x=combined_data[,c(gene_col_tumor)],base=2)<0] <- NA
      
      #on normal samples
    #  gene_col_normal <- paste(names(GENE_AND_ENSEMBL_ID)[GENE_AND_ENSEMBL_ID==GENE],"_",GENE,"_FPKM_normal",sep="")  
    #  combined_data[,c(gene_col_normal)][log(x=combined_data[,c(gene_col_normal)],base=2)<0] <- NA
      
      #divide non-NA by tumor purity
      #on tumor samples
  #    new_gene_col_tumor <- paste(names(GENE_AND_ENSEMBL_ID)[GENE_AND_ENSEMBL_ID==GENE],"_tumor_FPKM_adj",sep='')
   #   combined_data[,c(new_gene_col_tumor)] <- log(x=as.numeric(as.numeric(combined_data[,c(gene_col_tumor)]))/combined_data$Tumor_Purity,base=2)
      #on normal samples
  #    new_gene_col_normal <- paste(names(GENE_AND_ENSEMBL_ID)[GENE_AND_ENSEMBL_ID==GENE],"_normal_FPKM_adj",sep='')
   #   combined_data[,c(new_gene_col_normal)] <- log(x=as.numeric(as.numeric(combined_data[,c(gene_col_normal)]))/combined_data$Tumor_Purity,base=2)
  #  }
    
    #save dataset
    system(paste0('mkdir /home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/FPKM/',names(GENE_AND_ENSEMBL_ID)[GENE_AND_ENSEMBL_ID==GENE]))
    write.table(x=combined_data,file=paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/FPKM/',names(GENE_AND_ENSEMBL_ID)[GENE_AND_ENSEMBL_ID==GENE],'/',CANCER_NAME,'_',names(GENE_AND_ENSEMBL_ID)[GENE_AND_ENSEMBL_ID==GENE],'_FPKM_etc.txt'),sep='\t',row.names=FALSE)
}    
  
  
  