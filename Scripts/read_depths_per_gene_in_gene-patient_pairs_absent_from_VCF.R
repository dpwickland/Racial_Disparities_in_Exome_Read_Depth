#usage: Rscript --vanilla TCGA_gene_depths_pairs_missing_from_VCF.R BRCA TP53

args<-commandArgs(TRUE)
.libPaths()
# test whether both arguments present
if (length(args)!=2) {
  stop("Usage example: Rscript --vanilla TCGA_gene_depths_pairs_missing_from_VCF.R BRCA TP53", call.=FALSE)
} else if (length(args)==2) {
  # default output file
  CANCER_NAME <- toString(args[1])
  GENE <- toString(args[2])
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

#load overall depths data
combined_data_ID_and_bams_list <- read.table(paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/overall_depths/',CANCER_NAME,'_overall_depths_etc.txt'),header=TRUE)[,c(1,3,4,5)]
combined_data_ID_and_bams_list <- subset(combined_data_ID_and_bams_list, (read_depth_exome != "NA"))

#load gene depths data
cdriver_tumor_mutations <- read.table(paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/gene_depths_pairs_from_VCF/',GENE,'/',CANCER_NAME,'_',GENE,'_depths_etc.txt'),header=TRUE)[,c(1:5)]

#load all gene mutations data
#cdriver_tumor_mutations <- cdriver[,c(1:5)]
  #read.table(paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/gene_depths_pairs_from_VCF/',GENE,'/',CANCER_NAME,'_',GENE,'_all_detected_mutations.txt'),sep='\t')


############################################
#INDIVIDUAL GENE ANALYSIS: ALL PAIRWISE COMPARISONS BETWEEN EACH INDIVIDUAL'S BAM AND THE LIST OF ALL DETECTED SOMATIC VARIANTS IN ALL SAMPLES#
############################################

#1. Get list of just the unique loci; use these to query *ALL* individuals' bams at these coordinates, then order based on position
  cdriver_tumor_mutations_list <- unique(cdriver_tumor_mutations[,c(2,3,4,5)])
  cdriver_tumor_mutations_list <- cdriver_tumor_mutations_list[order(cdriver_tumor_mutations_list$mutation_position),]

#2. Create dataframe for all pairwise comparisons of patient BAM and gene variants detected in somatic VCFs
  #each case ID must be repeated same number of times as number of mutations
  #for example, if 496 individuals and 31 loci, then repeat 496 31 times
  #individual 1 will be repeated 31 times, then individual 2 will be repeated 31 times, etc.
  combined_data_ID_and_bams_list_repeated <- combined_data_ID_and_bams_list[rep(seq_len(nrow(combined_data_ID_and_bams_list)), each = nrow(cdriver_tumor_mutations_list)), ]
  
  #mutations list must be repeated for each case ID; repeated in same order as in original
  #ordered mutations list (based on position) of 31 mutations will be repeated 496 times
  cdriver_tumor_mutations_list_repeated <- cdriver_tumor_mutations_list[rep(1:nrow(cdriver_tumor_mutations_list),nrow(combined_data_ID_and_bams_list)),]
  
  #clean up environment
  rm(combined_data_ID_and_bams_list)
  rm(cdriver_tumor_mutations_list)
  
  #put it all together
  all_bams_all_mutations <- cbind(combined_data_ID_and_bams_list_repeated, cdriver_tumor_mutations_list_repeated)
  rm(combined_data_ID_and_bams_list_repeated)
  rm(cdriver_tumor_mutations_list_repeated)

#3. Read depth calculations: Count number of mapped reads of each all_bams_all_mutations mutated position in each individual
  variant_read_depth_bam_MPILEUP_col <- ncol(all_bams_all_mutations)+1
  ref_base_count_bam_col <- ncol(all_bams_all_mutations)+2
  alt_base_count_bam_col <- ncol(all_bams_all_mutations)+3
  #avg_read_mapping_quality_ref_bam_col <- ncol(all_bams_all_mutations)+4
  #avg_read_mapping_quality_alt_bam_col <- ncol(all_bams_all_mutations)+5
  avg_base_quality_ref_bam_col <- ncol(all_bams_all_mutations)+4
  avg_base_quality_alt_bam_col <- ncol(all_bams_all_mutations)+5
  
  
  for (i in 1:nrow(all_bams_all_mutations)){ 
    print(paste(i,"/",nrow(all_bams_all_mutations),sep=''))
    
    #calculate read depth at position using samtools view 
    #variant_read_depth_bam_VIEW <- system(paste("if [ -f ",sub(".bam",".bai",all_bams_all_mutations$bam_full_path_exome[i])," ]; then /research/bsi/tools/biotools/samtools/1.9/miniconda/bin/samtools view -c -F4 -F1024 -q 20",all_bams_all_mutations$bam_full_path_exome[i],paste(all_bams_all_mutations$chromosome[i],":",all_bams_all_mutations$mutation_position[i],"-",all_bams_all_mutations$mutation_position[i],"; else echo NA; fi",sep='')),intern=TRUE,ignore.stderr=TRUE)
    #all_bams_all_mutations[i,variant_read_depth_bam_VIEW_col] <- variant_read_depth_bam_VIEW
    
    #calculate read depth for REF, ALT and total using samtools mpileup
    #if (length(variant_read_depth_bam_VIEW!=0)){
    pileup <- system(paste("if [ -f ",sub(".bam",".bai",all_bams_all_mutations$bam_full_path_exome[i])," ]; then /research/bsi/tools/biotools/samtools/1.9/miniconda/bin/samtools mpileup ",all_bams_all_mutations$bam_full_path_exome[i],"--min-MQ 20 --min-BQ 0 --count-orphans --output-MQ --region ",paste(all_bams_all_mutations$chromosome[i],":",all_bams_all_mutations$mutation_position[i],"-",all_bams_all_mutations$mutation_position[i],"; else echo NA; fi",sep='')),intern=TRUE,ignore.stderr=TRUE)
    
    if (length(pileup)!=0 && pileup!= 'NA'){
      variant_read_depth_bam_MPILEUP <- strsplit(pileup,'\t')[[1]][4] 
      bases <- strsplit(pileup,'\t')[[1]][5]
      bases <- gsub("[^tcgaTCGA]+","",bases) #remove $ character from list of bases; this simply; don't need to do this for the qualities b/c no. qual values is already equal to no. bases
      
      ref_base_count_bam <- str_count(toupper(bases),toupper(as.character(all_bams_all_mutations[i,]$ref_allele_VCF))) #count number of reads with ref base in the pileup (convert all to uppercase)
      alt_base_count_bam <- str_count(toupper(bases),toupper(as.character(all_bams_all_mutations[i,]$alt_allele_VCF))) #count number of reads with alt base in the pileup (convert all to uppercase)
      
      all_bams_all_mutations[i,variant_read_depth_bam_MPILEUP_col] <- variant_read_depth_bam_MPILEUP
      all_bams_all_mutations[i,ref_base_count_bam_col] <- ref_base_count_bam
      all_bams_all_mutations[i,alt_base_count_bam_col] <- alt_base_count_bam
      
      ref_indices <- unlist(str_locate_all(toupper(bases),toupper(as.character(all_bams_all_mutations[i,]$ref_allele_VCF))))[1:ref_base_count_bam] #identify index positions of REF allele in the pileup; then get their corresponding quality scores; unlisting must be chopped in half b/c entire vector repeated once
      alt_indices <- unlist(str_locate_all(toupper(bases),toupper(as.character(all_bams_all_mutations[i,]$alt_allele_VCF))))[1:alt_base_count_bam] #identify index positions of REF allele in the pileup; then get their corresponding quality scores; unlisting must be chopped in half b/c entire vector repeated once
      
    #  read_mapping_quality <- asc(strsplit(pileup,'\t')[[1]][7])-33
    #  avg_read_mapping_quality_ref_bam <- round(mean(read_mapping_quality[ref_indices]),digits=2)
    #  avg_read_mapping_quality_alt_bam <- round(mean(read_mapping_quality[alt_indices]),digits=2) 
      
      base_quality <- asc(strsplit(pileup,'\t')[[1]][6])-33
      avg_base_quality_ref_bam <- round(mean(base_quality[ref_indices]),digits=2)
      avg_base_quality_alt_bam <- round(mean(base_quality[alt_indices]),digits=2)
      
      
      
    #  all_bams_all_mutations[i,avg_read_mapping_quality_ref_bam_col] <- avg_read_mapping_quality_ref_bam
    #  all_bams_all_mutations[i,avg_read_mapping_quality_alt_bam_col] <- avg_read_mapping_quality_alt_bam
      
      all_bams_all_mutations[i,avg_base_quality_ref_bam_col] <- avg_base_quality_ref_bam
      all_bams_all_mutations[i,avg_base_quality_alt_bam_col] <- avg_base_quality_alt_bam
      
    }
    
    else if (length(pileup)==0 || pileup=='NA'){
      all_bams_all_mutations[i,variant_read_depth_bam_MPILEUP_col] <- NA
      all_bams_all_mutations[i,ref_base_count_bam_col] <- NA
      all_bams_all_mutations[i,alt_base_count_bam_col] <- NA
      
      
  #    all_bams_all_mutations[i,avg_read_mapping_quality_ref_bam_col] <- NA
  #    all_bams_all_mutations[i,avg_read_mapping_quality_alt_bam_col] <- NA
      all_bams_all_mutations[i,avg_base_quality_ref_bam_col] <- NA
      all_bams_all_mutations[i,avg_base_quality_alt_bam_col] <- NA
      
    }
    
    #total_read_depth <- system(paste("/research/bsi/tools/biotools/samtools/1.9/miniconda/bin/samtools view -c -F4",all_bams_all_mutations$file_name_exome_bams[i]),intern=TRUE)
    #all_bams_all_mutations[i,total_read_depth_col] <- total_read_depth
    i <- i + 1 #increase row
  }
  
  colnames(all_bams_all_mutations)[variant_read_depth_bam_MPILEUP_col] <- 'variant_read_depth_bam_MPILEUP'
  colnames(all_bams_all_mutations)[ref_base_count_bam_col] <- 'ref_base_count_bam'
  colnames(all_bams_all_mutations)[alt_base_count_bam_col] <- 'alt_base_count_bam'
#  colnames(all_bams_all_mutations)[avg_read_mapping_quality_ref_bam_col] <- 'avg_read_mapping_quality_ref_bam'
#  colnames(all_bams_all_mutations)[avg_read_mapping_quality_alt_bam_col] <- 'avg_read_mapping_quality_alt_bam'
  colnames(all_bams_all_mutations)[avg_base_quality_ref_bam_col] <- 'avg_ref_base_quality_bam'
  colnames(all_bams_all_mutations)[avg_base_quality_alt_bam_col] <- 'avg_alt_base_quality_bam'


#4. Add additional column: alt allele concentration 
  all_bams_all_mutations$alt_allele_concent_bam <- as.numeric(all_bams_all_mutations$alt_base_count_bam) / as.numeric(all_bams_all_mutations$variant_read_depth_bam_MPILEUP)    
  
system(paste0("mkdir /home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/gene_depths_pairs_missing_from_VCF/",GENE))  

write.table(x=all_bams_all_mutations,file=paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/gene_depths_pairs_missing_from_VCF/',GENE,'/',CANCER_NAME,'_',GENE,'_depths_all_mutations_all_samples_etc.txt'),sep='\t',row.names=FALSE)  







