#usage: Rscript --vanilla TCGA_gene_depths_pairs_missing_from_VCF.R BRCA TP53

args<-commandArgs(TRUE)
.libPaths()
# test whether both arguments present
if (length(args)!=2) {
  stop("Usage example: Rscript --vanilla /home/mayo/m187735/s212975.Wickland_Immunomics/GitHub/Racial_Disparities_in_Exome_Read_Depth/Scripts/read_depths_per_gene_in_gene-patient_pairs_absent_from_MAF.R BRCA PIK3CA", call.=FALSE)
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

#load overall depths data, subsetting only Black and White samples
combined_data_ID_and_bams_list <- read.table(paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/overall_depths/',CANCER_NAME,'_overall_depths_etc.txt'),header=TRUE)[,c(1,3,4,5)]
combined_data_ID_and_bams_list <- subset(combined_data_ID_and_bams_list, (read_depth_exome != "NA" & (race == 'Black' | race =='White')))
combined_data_ID_and_bams_list <- droplevels(combined_data_ID_and_bams_list)

#load gene depths data
cdriver_tumor_mutations <- read.table(paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/gene_depths_pairs_from_MAF/',GENE,'/',CANCER_NAME,'_',GENE,'_depths_etc.txt'),header=TRUE)[,c(1,2,3,5,6)]

############################################
#INDIVIDUAL GENE ANALYSIS: ALL PAIRWISE COMPARISONS BETWEEN EACH INDIVIDUAL'S BAM AND THE LIST OF ALL DETECTED SOMATIC VARIANTS IN ALL SAMPLES#
############################################

#1. Get list of just the unique loci+allele; use these to query *ALL* individuals' bams at these coordinates, then order based on position
  cdriver_tumor_mutations_list <- unique(cdriver_tumor_mutations[,c(2,3,4,5)])
  cdriver_tumor_mutations_list <- cdriver_tumor_mutations_list[order(cdriver_tumor_mutations_list$Start_Position),]

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
  
  
  #remove those patient + allele combos that are present in the MAF file
  all_bams_all_mutations_absent_from_MAF <- all_bams_all_mutations[!interaction(all_bams_all_mutations[c('submitter_id','Chromosome','Start_Position','Tumor_Seq_Ref','Tumor_Seq_Alt')]) %in% interaction(cdriver_tumor_mutations[c('submitter_id','Chromosome','Start_Position','Tumor_Seq_Ref','Tumor_Seq_Alt')]),]
  
  

#3. Read depth calculations: Count number of mapped reads of each all_bams_all_mutations_absent_from_MAF mutated position in each individual
  variant_read_depth_bam_MPILEUP_col <- ncol(all_bams_all_mutations_absent_from_MAF)+1
  ref_base_count_bam_col <- ncol(all_bams_all_mutations_absent_from_MAF)+2
  alt_base_count_bam_col <- ncol(all_bams_all_mutations_absent_from_MAF)+3
  #avg_read_mapping_quality_ref_bam_col <- ncol(all_bams_all_mutations_absent_from_MAF)+4
  #avg_read_mapping_quality_alt_bam_col <- ncol(all_bams_all_mutations_absent_from_MAF)+5
  avg_base_quality_ref_bam_col <- ncol(all_bams_all_mutations_absent_from_MAF)+4
  avg_base_quality_alt_bam_col <- ncol(all_bams_all_mutations_absent_from_MAF)+5
  
  
  for (i in 1:nrow(all_bams_all_mutations_absent_from_MAF)){ 
    print(paste(i,"/",nrow(all_bams_all_mutations_absent_from_MAF),sep=''))
    

    #calculate read depth for REF, ALT and total using samtools mpileup
    pileup <- system(paste("if [ -f ",sub(".bam",".bai",all_bams_all_mutations_absent_from_MAF$bam_full_path_exome[i])," ]; then /research/bsi/tools/biotools/samtools/1.9/miniconda/bin/samtools mpileup ",all_bams_all_mutations_absent_from_MAF$bam_full_path_exome[i],"--min-MQ 20 --min-BQ 0 --count-orphans --output-MQ --region ",paste(all_bams_all_mutations_absent_from_MAF$Chromosome[i],":",all_bams_all_mutations_absent_from_MAF$Start_Position[i],"-",all_bams_all_mutations_absent_from_MAF$Start_Position[i],"; else echo NA; fi",sep='')),intern=TRUE,ignore.stderr=TRUE)
    
    if (length(pileup)!=0 && pileup!= 'NA'){
      variant_read_depth_bam_MPILEUP <- strsplit(pileup,'\t')[[1]][4] 
      bases <- strsplit(pileup,'\t')[[1]][5]
      bases <- gsub("[^tcgaTCGA]+","",bases) #remove $ character from list of bases; this simply; don't need to do this for the qualities b/c no. qual values is already equal to no. bases
      
      ref_base_count_bam <- str_count(toupper(bases),toupper(as.character(all_bams_all_mutations_absent_from_MAF[i,]$Tumor_Seq_Ref))) #count number of reads with ref base in the pileup (convert all to uppercase)
      alt_base_count_bam <- str_count(toupper(bases),toupper(as.character(all_bams_all_mutations_absent_from_MAF[i,]$Tumor_Seq_Alt))) #count number of reads with alt base in the pileup (convert all to uppercase)
      
      all_bams_all_mutations_absent_from_MAF[i,variant_read_depth_bam_MPILEUP_col] <- variant_read_depth_bam_MPILEUP
      all_bams_all_mutations_absent_from_MAF[i,ref_base_count_bam_col] <- ref_base_count_bam
      all_bams_all_mutations_absent_from_MAF[i,alt_base_count_bam_col] <- alt_base_count_bam
      
      ref_indices <- unlist(str_locate_all(toupper(bases),toupper(as.character(all_bams_all_mutations_absent_from_MAF[i,]$Tumor_Seq_Ref))))[1:ref_base_count_bam] #identify index positions of REF allele in the pileup; then get their corresponding quality scores; unlisting must be chopped in half b/c entire vector repeated once
      alt_indices <- unlist(str_locate_all(toupper(bases),toupper(as.character(all_bams_all_mutations_absent_from_MAF[i,]$Tumor_Seq_Alt))))[1:alt_base_count_bam] #identify index positions of REF allele in the pileup; then get their corresponding quality scores; unlisting must be chopped in half b/c entire vector repeated once
      
    #  read_mapping_quality <- asc(strsplit(pileup,'\t')[[1]][7])-33
    #  avg_read_mapping_quality_ref_bam <- round(mean(read_mapping_quality[ref_indices]),digits=2)
    #  avg_read_mapping_quality_alt_bam <- round(mean(read_mapping_quality[alt_indices]),digits=2) 
      
      base_quality <- asc(strsplit(pileup,'\t')[[1]][6])-33
      avg_base_quality_ref_bam <- round(mean(base_quality[ref_indices]),digits=2)
      avg_base_quality_alt_bam <- round(mean(base_quality[alt_indices]),digits=2)
      
      
      
    #  all_bams_all_mutations_absent_from_MAF[i,avg_read_mapping_quality_ref_bam_col] <- avg_read_mapping_quality_ref_bam
    #  all_bams_all_mutations_absent_from_MAF[i,avg_read_mapping_quality_alt_bam_col] <- avg_read_mapping_quality_alt_bam
      
      all_bams_all_mutations_absent_from_MAF[i,avg_base_quality_ref_bam_col] <- avg_base_quality_ref_bam
      all_bams_all_mutations_absent_from_MAF[i,avg_base_quality_alt_bam_col] <- avg_base_quality_alt_bam
      
    }
    
    else if (length(pileup)==0 || pileup=='NA'){
      all_bams_all_mutations_absent_from_MAF[i,variant_read_depth_bam_MPILEUP_col] <- NA
      all_bams_all_mutations_absent_from_MAF[i,ref_base_count_bam_col] <- NA
      all_bams_all_mutations_absent_from_MAF[i,alt_base_count_bam_col] <- NA
      
      
  #    all_bams_all_mutations_absent_from_MAF[i,avg_read_mapping_quality_ref_bam_col] <- NA
  #    all_bams_all_mutations_absent_from_MAF[i,avg_read_mapping_quality_alt_bam_col] <- NA
      all_bams_all_mutations_absent_from_MAF[i,avg_base_quality_ref_bam_col] <- NA
      all_bams_all_mutations_absent_from_MAF[i,avg_base_quality_alt_bam_col] <- NA
      
    }
    
    #total_read_depth <- system(paste("/research/bsi/tools/biotools/samtools/1.9/miniconda/bin/samtools view -c -F4",all_bams_all_mutations_absent_from_MAF$file_name_exome_bams[i]),intern=TRUE)
    #all_bams_all_mutations_absent_from_MAF[i,total_read_depth_col] <- total_read_depth
    i <- i + 1 #increase row
  }
  
  colnames(all_bams_all_mutations_absent_from_MAF)[variant_read_depth_bam_MPILEUP_col] <- 'variant_read_depth_bam_MPILEUP'
  colnames(all_bams_all_mutations_absent_from_MAF)[ref_base_count_bam_col] <- 'ref_base_count_bam'
  colnames(all_bams_all_mutations_absent_from_MAF)[alt_base_count_bam_col] <- 'alt_base_count_bam'
  #  colnames(all_bams_all_mutations_absent_from_MAF)[avg_read_mapping_quality_ref_bam_col] <- 'avg_read_mapping_quality_ref_bam'
  #  colnames(all_bams_all_mutations_absent_from_MAF)[avg_read_mapping_quality_alt_bam_col] <- 'avg_read_mapping_quality_alt_bam'
  colnames(all_bams_all_mutations_absent_from_MAF)[avg_base_quality_ref_bam_col] <- 'avg_ref_base_quality_bam'
  colnames(all_bams_all_mutations_absent_from_MAF)[avg_base_quality_alt_bam_col] <- 'avg_alt_base_quality_bam'


#4. Add additional column: alt allele concentration 
  all_bams_all_mutations_absent_from_MAF$alt_allele_concent_bam <- as.numeric(all_bams_all_mutations_absent_from_MAF$alt_base_count_bam) / as.numeric(all_bams_all_mutations_absent_from_MAF$variant_read_depth_bam_MPILEUP)    
  
system(paste0("mkdir -p /home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/gene_depths_pairs_missing_from_MAF/",GENE))  

write.table(x=all_bams_all_mutations_absent_from_MAF,file=paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/gene_depths_pairs_missing_from_MAF/',GENE,'/',CANCER_NAME,'_',GENE,'_depths_all_mutations_all_samples_absent_from_MAF_etc.txt'),sep='\t',row.names=FALSE)  






