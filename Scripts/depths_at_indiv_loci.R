GENE <- 'PIK3CA'


CANCER_LIST <- c('BRCA','LUAD','UCEC','PRAD','COAD')

all_cases_all_cancers <- data.frame()


for (CANCER_NAME in CANCER_LIST){
  #############################################################
  #LOAD DATA
  #############################################################
  
  #load gene depths data
  cases_with_MAF_mutations <- read.table(paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/gene_depths_pairs_from_MAF/',GENE,'/',CANCER_NAME,'_',GENE,'_depths_etc.txt'),header=TRUE)
 cases_with_MAF_mutations$Cancer <- CANCER_NAME
 cases_with_MAF_mutations$Gene <- GENE
   cases_with_MAF_mutations$`Patient+SNP in MAF?` <- "Yes"
   cases_with_MAF_mutations <- subset(cases_with_MAF_mutations, (race=='Black' | race=='White'))
  
   cases_with_MAF_mutations <- cases_with_MAF_mutations[,c(1,28,2,3,5,6,15,16,17,23,24,25,32,38,39,40)]
   names(cases_with_MAF_mutations)[14] <- 'SNP in at least 1 case in MAF?'

   colnames(cases_with_MAF_mutations)[7] <- "MAF_depth"
   names(cases_with_MAF_mutations)[8] <- "MAF_ref_depth"
   names(cases_with_MAF_mutations)[9] <- "MAF_alt_depth"
   colnames(cases_with_MAF_mutations)[10] <- "BAM_depth"
   names(cases_with_MAF_mutations)[11] <- "BAM_ref_depth"
   names(cases_with_MAF_mutations)[12] <- "BAM_alt_depth"
   cases_with_MAF_mutations <- subset(cases_with_MAF_mutations, (submitter_id!='NA'))
   
  #load all bams all mutations in gene data
  cases_without_MAF_mutations <- read.table(paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/gene_depths_pairs_missing_from_MAF/',GENE,'/',CANCER_NAME,'_',GENE,'_depths_all_mutations_all_samples_absent_from_MAF_etc.txt'),header=TRUE)
  cases_without_MAF_mutations$Cancer <- CANCER_NAME
  cases_without_MAF_mutations$Gene <- GENE
  cases_without_MAF_mutations$`Patient+SNP in MAF?` <- "No"
  
  cases_without_MAF_mutations <- cases_without_MAF_mutations[,c(1,2,5,6,7,8,9,10,11,15,16,17)]
  names(cases_without_MAF_mutations)[7] <- "BAM_depth"
  names(cases_without_MAF_mutations)[8] <- "BAM_ref_depth"
  names(cases_without_MAF_mutations)[9] <- "BAM_alt_depth"
  
  cases_without_MAF_mutations$MAF_depth <- NA
  cases_without_MAF_mutations$MAF_ref_depth <- NA
  cases_without_MAF_mutations$MAF_alt_depth <- NA
  
  cases_without_MAF_mutations_but_SNPs_in_MAF <- cases_without_MAF_mutations[paste0(cases_without_MAF_mutations$Start_Position, cases_without_MAF_mutations$Tumor_Seq_Alt) %in% paste0(cases_with_MAF_mutations$Start_Position, cases_with_MAF_mutations$Tumor_Seq_Alt),]
  cases_without_MAF_mutations_but_SNPs_in_MAF$`SNP in at least 1 case in MAF?` <- "Yes"
  
  `%notin%` = Negate(`%in%`)
  cases_without_MAF_mutations_and_SNPs_not_in_MAF <- cases_without_MAF_mutations[paste0(cases_without_MAF_mutations$Start_Position, cases_without_MAF_mutations$Tumor_Seq_Alt) %notin% paste0(cases_with_MAF_mutations$Start_Position, cases_with_MAF_mutations$Tumor_Seq_Alt),]
  cases_without_MAF_mutations_and_SNPs_not_in_MAF$`SNP in at least 1 case in MAF?` <- "No"

  
  all_cases <- rbind(cases_with_MAF_mutations, cases_without_MAF_mutations_but_SNPs_in_MAF, cases_without_MAF_mutations_and_SNPs_not_in_MAF)
  
  colnames(all_cases)[4] <- "ClinVar Position"
  
  #load VCF depths data
  cases_with_VCF_mutations <- read.table(paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/gene_depths_pairs_from_VCF/',GENE,'/',CANCER_NAME,'_',GENE,'_depths_etc.txt'),header=TRUE)
  cases_with_VCF_mutations <- cases_with_VCF_mutations[,c('submitter_id','chromosome','mutation_position','ref_allele_VCF','alt_allele_VCF','AD_total_VCF','AD_ref_VCF','AD_alt_VCF')]
  
  colnames(cases_with_VCF_mutations) <- c('submitter_id','Chromosome','ClinVar Position','Tumor_Seq_Ref','Tumor_Seq_Alt','VCF_depth','VCF_ref_depth','VCF_alt_depth')
  
  all_cases <- merge(all_cases, cases_with_VCF_mutations, by=c('submitter_id','Chromosome','ClinVar Position','Tumor_Seq_Ref','Tumor_Seq_Alt'), all.x=TRUE)
  
  
  all_cases_all_cancers <- rbind(all_cases_all_cancers, all_cases)
  

}

#order by position
all_cases_all_cancers <- all_cases_all_cancers[order(all_cases_all_cancers$`ClinVar Position`,all_cases_all_cancers$race),]

#reorder columns
all_cases_all_cancers <- all_cases_all_cancers[,c(1,6,2:5,7:9,17:19,10:16)]

write.table(x=all_cases_all_cancers, file=paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/read_depth_table_',GENE,'.txt'),sep='\t',row.names=FALSE)  


#table for number of mutations, across cancers, called in B + W
race_chrom_position_alt <- all_cases_all_cancers[, c("race","Chromosome","Cancer","ClinVar Position","Tumor_Seq_Alt", "SNP in at least 1 case in MAF?","Patient+SNP in MAF?")]

positions_called_in_at_least_one_case_across_all  <- subset(race_chrom_position_alt, (`SNP in at least 1 case in MAF?`=='Yes'))$`ClinVar Position`
race_chrom_position_alt <- race_chrom_position_alt[race_chrom_position_alt$`ClinVar Position` %in% positions_called_in_at_least_one_case_across_all,]

summary_table <- race_chrom_position_alt %>% group_by(`ClinVar Position`,race,Tumor_Seq_Alt,`Patient+SNP in MAF?`) %>% summarise(n=n())
write.table(x=summary_table, file=paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/read_depth_summary_table_',GENE,'.txt'), sep='\t', row.names=FALSE)





                                                             