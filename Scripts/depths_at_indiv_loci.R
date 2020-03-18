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

   colnames(cases_with_MAF_mutations)[7] <- "SNP_depth_maf"
   names(cases_with_MAF_mutations)[8] <- "SNP_ref_depth_maf"
   names(cases_with_MAF_mutations)[9] <- "SNP_alt_depth_maf"
   colnames(cases_with_MAF_mutations)[10] <- "SNP_depth_bam"
   names(cases_with_MAF_mutations)[11] <- "SNP_ref_depth_bam"
   names(cases_with_MAF_mutations)[12] <- "SNP_alt_depth_bam"
   cases_with_MAF_mutations <- subset(cases_with_MAF_mutations, (submitter_id!='NA'))
   
  #load all bams all mutations in gene data
  cases_without_MAF_mutations <- read.table(paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/gene_depths_pairs_missing_from_MAF/',GENE,'/',CANCER_NAME,'_',GENE,'_depths_all_mutations_all_samples_absent_from_MAF_etc.txt'),header=TRUE)
  cases_without_MAF_mutations$Cancer <- CANCER_NAME
  cases_without_MAF_mutations$Gene <- GENE
  cases_without_MAF_mutations$`Patient+SNP in MAF?` <- "No"
  
  cases_without_MAF_mutations <- cases_without_MAF_mutations[,c(1,2,5,6,7,8,9,10,11,15,16,17)]
  names(cases_without_MAF_mutations)[7] <- "SNP_depth_bam"
  names(cases_without_MAF_mutations)[8] <- "SNP_ref_depth_bam"
  names(cases_without_MAF_mutations)[9] <- "SNP_alt_depth_bam"
  
  cases_without_MAF_mutations$SNP_depth_maf <- NA
  cases_without_MAF_mutations$SNP_ref_depth_maf <- NA
  cases_without_MAF_mutations$SNP_alt_depth_maf <- NA
  
  cases_without_MAF_mutations_but_SNPs_in_MAF <- cases_without_MAF_mutations[paste0(cases_without_MAF_mutations$Start_Position, cases_without_MAF_mutations$Tumor_Seq_Alt) %in% paste0(cases_with_MAF_mutations$Start_Position, cases_with_MAF_mutations$Tumor_Seq_Alt),]
  cases_without_MAF_mutations_but_SNPs_in_MAF$`SNP in at least 1 case in MAF?` <- "Yes"
  
  `%notin%` = Negate(`%in%`)
  cases_without_MAF_mutations_and_SNPs_not_in_MAF <- cases_without_MAF_mutations[paste0(cases_without_MAF_mutations$Start_Position, cases_without_MAF_mutations$Tumor_Seq_Alt) %notin% paste0(cases_with_MAF_mutations$Start_Position, cases_with_MAF_mutations$Tumor_Seq_Alt),]
  cases_without_MAF_mutations_and_SNPs_not_in_MAF$`SNP in at least 1 case in MAF?` <- "No"

  
  all_cases <- rbind(cases_with_MAF_mutations, cases_without_MAF_mutations_but_SNPs_in_MAF, cases_without_MAF_mutations_and_SNPs_not_in_MAF)
  
  colnames(all_cases)[4] <- "ClinVar Position"
  
  all_cases_all_cancers <- rbind(all_cases_all_cancers, all_cases)
  

}

#order by position
all_cases_all_cancers <- all_cases_all_cancers[order(all_cases_all_cancers$`ClinVar Position`,all_cases_all_cancers$race),]

write.table(x=all_cases_all_cancers, file=paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/read_depth_table_',GENE,'.txt'),sep='\t',row.names=FALSE)  


#table for number of mutations, across cancers, called in B + W
race_chrom_position_alt <- all_cases_all_cancers[, c("race","Chromosome","Cancer","ClinVar Position","Tumor_Seq_Alt", "SNP in at least 1 case in MAF?","Patient+SNP in MAF?")]

positions_called_in_at_least_one_case_across_all  <- subset(race_chrom_position_alt, (`SNP in at least 1 case in MAF?`=='Yes'))$`ClinVar Position`
race_chrom_position_alt <- race_chrom_position_alt[race_chrom_position_alt$`ClinVar Position` %in% positions_called_in_at_least_one_case_across_all,]

summary_table <- race_chrom_position_alt %>% group_by(`ClinVar Position`,race,Tumor_Seq_Alt,`Patient+SNP in MAF?`) %>% summarise(n=n())
write.table(x=summary_table, file=paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/read_depth_summary_table_',GENE,'.txt'), sep='\t', row.names=FALSE)





                                                             