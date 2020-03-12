  
  CANCER_LIST <- c('BRCA','LUAD','UCEC','PRAD','COAD')
  
  GENE <- 'TP53'
  
  Black_White_cases_with_MAF_mutations <- data.frame()
  Black_White_cases_without_MAF_mutations <- data.frame()
  
  
  for (CANCER_NAME in CANCER_LIST){
    #############################################################
    #LOAD DATA
    #############################################################
    
    #load gene depths data
    cases_with_MAF_mutations <- read.table(paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/gene_depths_pairs_from_MAF/',GENE,'/',CANCER_NAME,'_',GENE,'_depths_etc.txt'),header=TRUE)
    
    #load all bams all mutations in gene data
    cases_without_MAF_mutations <- read.table(paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/gene_depths_pairs_missing_from_MAF/',GENE,'/',CANCER_NAME,'_',GENE,'_depths_all_mutations_all_samples_absent_from_MAF_etc.txt'),header=TRUE)
    cases_without_MAF_mutations$Cancer <- CANCER_NAME
  
    #order by race
    cases_with_MAF_mutations <- subset(cases_with_MAF_mutations, (race == 'White' | race == 'Black'))
    cases_with_MAF_mutations$race <- factor(cases_with_MAF_mutations$race, levels = c("White","Black")) 
    cases_without_MAF_mutations$race <- factor(cases_without_MAF_mutations$race, levels = c("White","Black")) 
    
    #create separate dataframes for white and black
    black_cases_with_MAF_mutations <- subset(cases_with_MAF_mutations, (race=='Black'))
    black_cases_without_MAF_mutations <- subset(cases_without_MAF_mutations, (race=='Black'))
    
    white_cases_with_MAF_mutations <- subset(cases_with_MAF_mutations, (race=='White'))
    white_cases_without_MAF_mutations <- subset(cases_without_MAF_mutations, (race=='White'))
    
    #load combined data to determine number of samples for each cancer
    combined_data <- read.table(paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/overall_depths/',CANCER_NAME,'_overall_depths_etc.txt'),header=TRUE)
    
    if ("Black" %in% cases_with_MAF_mutations$race){
    
      rm(cases_with_MAF_mutations)
      rm(cases_without_MAF_mutations)
      
      #############################################################
      #REARRANGE DATA AND CALCULATE NUMBER/PERCENT OF PATIENTS WITH DIFFERING COVERAGES: CASES WITH MAF MUTATIONS
      #############################################################
      
      ###
      #BLACK CASES
      ###
      
        #BASED ON TOTAL DEPTH (REF + ALT)
        black_cases_with_MAF_mutations_reshaped_total_depth <- black_cases_with_MAF_mutations[,c(1,2,3,5,6,23,24,25,28,32)]
        black_cases_with_MAF_mutations_reshaped_total_depth <- dcast(black_cases_with_MAF_mutations_reshaped_total_depth, Cancer+race+Chromosome+Start_Position+Tumor_Seq_Ref+Tumor_Seq_Alt ~ submitter_id, value.var = 'variant_read_depth_bam_MPILEUP')
        
        if (ncol(black_cases_with_MAF_mutations_reshaped_total_depth) > 7){
        black_cases_with_MAF_mutations_reshaped_total_depth$total_depth_1x <- rowSums(na.rm=TRUE,black_cases_with_MAF_mutations_reshaped_total_depth[, 7:ncol(black_cases_with_MAF_mutations_reshaped_total_depth)] >=1) / nrow(subset(combined_data, (race=='Black'))) #length(7:ncol(black_cases_with_MAF_mutations_reshaped_total_depth))
        black_cases_with_MAF_mutations_reshaped_total_depth$total_depth_5x <- rowSums(na.rm=TRUE,black_cases_with_MAF_mutations_reshaped_total_depth[, 7:ncol(black_cases_with_MAF_mutations_reshaped_total_depth)] >=5) / nrow(subset(combined_data, (race=='Black')))
        black_cases_with_MAF_mutations_reshaped_total_depth$total_depth_10x <- rowSums(na.rm=TRUE,black_cases_with_MAF_mutations_reshaped_total_depth[, 7:ncol(black_cases_with_MAF_mutations_reshaped_total_depth)] >=10) / nrow(subset(combined_data, (race=='Black')))
        black_cases_with_MAF_mutations_reshaped_total_depth$total_depth_20x <- rowSums(na.rm=TRUE,black_cases_with_MAF_mutations_reshaped_total_depth[, 7:ncol(black_cases_with_MAF_mutations_reshaped_total_depth)] >=20) / nrow(subset(combined_data, (race=='Black')))
        black_cases_with_MAF_mutations_reshaped_total_depth$total_depth_50x <- rowSums(na.rm=TRUE,black_cases_with_MAF_mutations_reshaped_total_depth[, 7:ncol(black_cases_with_MAF_mutations_reshaped_total_depth)] >=50) / nrow(subset(combined_data, (race=='Black')))
        black_cases_with_MAF_mutations_reshaped_total_depth$total_depth_100x <- rowSums(na.rm=TRUE,black_cases_with_MAF_mutations_reshaped_total_depth[, 7:ncol(black_cases_with_MAF_mutations_reshaped_total_depth)] >=100) / nrow(subset(combined_data, (race=='Black')))
        
        black_cases_with_MAF_mutations_reshaped_total_depth <- black_cases_with_MAF_mutations_reshaped_total_depth[,c(1:6,(ncol(black_cases_with_MAF_mutations_reshaped_total_depth)-4):ncol(black_cases_with_MAF_mutations_reshaped_total_depth))]
        
        
        #BASED ON REF DEPTH
        black_cases_with_MAF_mutations_reshaped_ref_depth <- black_cases_with_MAF_mutations[,c(1,2,3,5,6,23,24,25,28,32)]
        black_cases_with_MAF_mutations_reshaped_ref_depth <- dcast(black_cases_with_MAF_mutations_reshaped_ref_depth, Cancer+race+Chromosome+Start_Position+Tumor_Seq_Ref+Tumor_Seq_Alt ~ submitter_id, value.var = 'tumor_ref_base_count_bam')
        black_cases_with_MAF_mutations_reshaped_ref_depth$ref_depth_1x <- rowSums(na.rm=TRUE,black_cases_with_MAF_mutations_reshaped_ref_depth[, 7:ncol(black_cases_with_MAF_mutations_reshaped_ref_depth)] >=1) / nrow(subset(combined_data, (race=='Black')))
        black_cases_with_MAF_mutations_reshaped_ref_depth$ref_depth_5x <- rowSums(na.rm=TRUE,black_cases_with_MAF_mutations_reshaped_ref_depth[, 7:ncol(black_cases_with_MAF_mutations_reshaped_ref_depth)] >=5) / nrow(subset(combined_data, (race=='Black')))
        black_cases_with_MAF_mutations_reshaped_ref_depth$ref_depth_10x <- rowSums(na.rm=TRUE,black_cases_with_MAF_mutations_reshaped_ref_depth[, 7:ncol(black_cases_with_MAF_mutations_reshaped_ref_depth)] >=10) / nrow(subset(combined_data, (race=='Black')))
        black_cases_with_MAF_mutations_reshaped_ref_depth$ref_depth_20x <- rowSums(na.rm=TRUE,black_cases_with_MAF_mutations_reshaped_ref_depth[, 7:ncol(black_cases_with_MAF_mutations_reshaped_ref_depth)] >=20) / nrow(subset(combined_data, (race=='Black')))
        black_cases_with_MAF_mutations_reshaped_ref_depth$ref_depth_50x <- rowSums(na.rm=TRUE,black_cases_with_MAF_mutations_reshaped_ref_depth[, 7:ncol(black_cases_with_MAF_mutations_reshaped_ref_depth)] >=50) / nrow(subset(combined_data, (race=='Black')))
        black_cases_with_MAF_mutations_reshaped_ref_depth$ref_depth_100x <- rowSums(na.rm=TRUE,black_cases_with_MAF_mutations_reshaped_ref_depth[, 7:ncol(black_cases_with_MAF_mutations_reshaped_ref_depth)] >=100) / nrow(subset(combined_data, (race=='Black')))
        
        black_cases_with_MAF_mutations_reshaped_ref_depth <- black_cases_with_MAF_mutations_reshaped_ref_depth[,c(1:6,(ncol(black_cases_with_MAF_mutations_reshaped_ref_depth)-4):ncol(black_cases_with_MAF_mutations_reshaped_ref_depth))]
        
        
        #BASED ON ALT DEPTH
        black_cases_with_MAF_mutations_reshaped_alt_depth <- black_cases_with_MAF_mutations[,c(1,2,3,5,6,23,24,25,28,32)]
        black_cases_with_MAF_mutations_reshaped_alt_depth <- dcast(black_cases_with_MAF_mutations_reshaped_alt_depth, Cancer+race+Chromosome+Start_Position+Tumor_Seq_Ref+Tumor_Seq_Alt ~ submitter_id, value.var = 'tumor_alt_base_count_bam')
        black_cases_with_MAF_mutations_reshaped_alt_depth$alt_depth_1x <- rowSums(na.rm=TRUE,black_cases_with_MAF_mutations_reshaped_alt_depth[, 7:ncol(black_cases_with_MAF_mutations_reshaped_alt_depth)] >=1) / nrow(subset(combined_data, (race=='Black')))
        black_cases_with_MAF_mutations_reshaped_alt_depth$alt_depth_5x <- rowSums(na.rm=TRUE,black_cases_with_MAF_mutations_reshaped_alt_depth[, 7:ncol(black_cases_with_MAF_mutations_reshaped_alt_depth)] >=5) / nrow(subset(combined_data, (race=='Black')))
        black_cases_with_MAF_mutations_reshaped_alt_depth$alt_depth_10x <- rowSums(na.rm=TRUE,black_cases_with_MAF_mutations_reshaped_alt_depth[, 7:ncol(black_cases_with_MAF_mutations_reshaped_alt_depth)] >=10) / nrow(subset(combined_data, (race=='Black')))
        black_cases_with_MAF_mutations_reshaped_alt_depth$alt_depth_20x <- rowSums(na.rm=TRUE,black_cases_with_MAF_mutations_reshaped_alt_depth[, 7:ncol(black_cases_with_MAF_mutations_reshaped_alt_depth)] >=20) / nrow(subset(combined_data, (race=='Black')))
        black_cases_with_MAF_mutations_reshaped_alt_depth$alt_depth_50x <- rowSums(na.rm=TRUE,black_cases_with_MAF_mutations_reshaped_alt_depth[, 7:ncol(black_cases_with_MAF_mutations_reshaped_alt_depth)] >=50) / nrow(subset(combined_data, (race=='Black')))
        black_cases_with_MAF_mutations_reshaped_alt_depth$alt_depth_100x <- rowSums(na.rm=TRUE,black_cases_with_MAF_mutations_reshaped_alt_depth[, 7:ncol(black_cases_with_MAF_mutations_reshaped_alt_depth)] >=100) / nrow(subset(combined_data, (race=='Black')))
        
        #MERGE DATA
        black_cases_with_MAF_mutations_reshaped_alt_depth <- black_cases_with_MAF_mutations_reshaped_alt_depth[,c(1:6,(ncol(black_cases_with_MAF_mutations_reshaped_alt_depth)-4):ncol(black_cases_with_MAF_mutations_reshaped_alt_depth))]
        black_cases_with_MAF_mutations_reshaped <- merge(black_cases_with_MAF_mutations_reshaped_total_depth, black_cases_with_MAF_mutations_reshaped_ref_depth, by=c('Cancer','race','Chromosome','Start_Position','Tumor_Seq_Ref','Tumor_Seq_Alt'))
          black_cases_with_MAF_mutations_reshaped <- merge(black_cases_with_MAF_mutations_reshaped, black_cases_with_MAF_mutations_reshaped_alt_depth, by=c('Cancer','race','Chromosome','Start_Position','Tumor_Seq_Ref','Tumor_Seq_Alt'))
        black_cases_with_MAF_mutations_reshaped <- black_cases_with_MAF_mutations_reshaped[order(black_cases_with_MAF_mutations_reshaped$Start_Position),]
        black_cases_with_MAF_mutations_reshaped[,c(7:ncol(black_cases_with_MAF_mutations_reshaped))] <- round(black_cases_with_MAF_mutations_reshaped[,c(7:ncol(black_cases_with_MAF_mutations_reshaped))]*100, digits=2)
        
        #ADD PERCENTAGE SIGN
        black_cases_with_MAF_mutations_reshaped[,c(7:ncol(black_cases_with_MAF_mutations_reshaped))] <- lapply(black_cases_with_MAF_mutations_reshaped[,c(7:ncol(black_cases_with_MAF_mutations_reshaped))], function(x) paste0(x,"%"))
        rm(black_cases_with_MAF_mutations_reshaped_alt_depth, black_cases_with_MAF_mutations_reshaped_ref_depth, black_cases_with_MAF_mutations_reshaped_total_depth)
      }
        ###
        #WHITE CASES
        ###
        
        #BASED ON TOTAL DEPTH (REF + ALT)
        white_cases_with_MAF_mutations_reshaped_total_depth <- white_cases_with_MAF_mutations[,c(1,2,3,5,6,23,24,25,28,32)]
        white_cases_with_MAF_mutations_reshaped_total_depth <- dcast(white_cases_with_MAF_mutations_reshaped_total_depth, Cancer+race+Chromosome+Start_Position+Tumor_Seq_Ref+Tumor_Seq_Alt ~ submitter_id, value.var = 'variant_read_depth_bam_MPILEUP')
        
        if (ncol(white_cases_with_MAF_mutations_reshaped_total_depth) > 7){
          white_cases_with_MAF_mutations_reshaped_total_depth$total_depth_1x <- rowSums(na.rm=TRUE,white_cases_with_MAF_mutations_reshaped_total_depth[, 7:ncol(white_cases_with_MAF_mutations_reshaped_total_depth)] >=1) / nrow(subset(combined_data, (race=='White')))
          white_cases_with_MAF_mutations_reshaped_total_depth$total_depth_5x <- rowSums(na.rm=TRUE,white_cases_with_MAF_mutations_reshaped_total_depth[, 7:ncol(white_cases_with_MAF_mutations_reshaped_total_depth)] >=5) / nrow(subset(combined_data, (race=='White')))
          white_cases_with_MAF_mutations_reshaped_total_depth$total_depth_10x <- rowSums(na.rm=TRUE,white_cases_with_MAF_mutations_reshaped_total_depth[, 7:ncol(white_cases_with_MAF_mutations_reshaped_total_depth)] >=10) / nrow(subset(combined_data, (race=='White')))
          white_cases_with_MAF_mutations_reshaped_total_depth$total_depth_20x <- rowSums(na.rm=TRUE,white_cases_with_MAF_mutations_reshaped_total_depth[, 7:ncol(white_cases_with_MAF_mutations_reshaped_total_depth)] >=20) / nrow(subset(combined_data, (race=='White')))
          white_cases_with_MAF_mutations_reshaped_total_depth$total_depth_50x <- rowSums(na.rm=TRUE,white_cases_with_MAF_mutations_reshaped_total_depth[, 7:ncol(white_cases_with_MAF_mutations_reshaped_total_depth)] >=50) / nrow(subset(combined_data, (race=='White')))
          white_cases_with_MAF_mutations_reshaped_total_depth$total_depth_100x <- rowSums(na.rm=TRUE,white_cases_with_MAF_mutations_reshaped_total_depth[, 7:ncol(white_cases_with_MAF_mutations_reshaped_total_depth)] >=100) / nrow(subset(combined_data, (race=='White')))
          
          white_cases_with_MAF_mutations_reshaped_total_depth <- white_cases_with_MAF_mutations_reshaped_total_depth[,c(1:6,(ncol(white_cases_with_MAF_mutations_reshaped_total_depth)-4):ncol(white_cases_with_MAF_mutations_reshaped_total_depth))]
          
          
          #BASED ON REF DEPTH
          white_cases_with_MAF_mutations_reshaped_ref_depth <- white_cases_with_MAF_mutations[,c(1,2,3,5,6,23,24,25,28,32)]
          white_cases_with_MAF_mutations_reshaped_ref_depth <- dcast(white_cases_with_MAF_mutations_reshaped_ref_depth, Cancer+race+Chromosome+Start_Position+Tumor_Seq_Ref+Tumor_Seq_Alt ~ submitter_id, value.var = 'tumor_ref_base_count_bam')
          white_cases_with_MAF_mutations_reshaped_ref_depth$ref_depth_1x <- rowSums(na.rm=TRUE,white_cases_with_MAF_mutations_reshaped_ref_depth[, 7:ncol(white_cases_with_MAF_mutations_reshaped_ref_depth)] >=1) / nrow(subset(combined_data, (race=='White')))
          white_cases_with_MAF_mutations_reshaped_ref_depth$ref_depth_5x <- rowSums(na.rm=TRUE,white_cases_with_MAF_mutations_reshaped_ref_depth[, 7:ncol(white_cases_with_MAF_mutations_reshaped_ref_depth)] >=5) / nrow(subset(combined_data, (race=='White')))
          white_cases_with_MAF_mutations_reshaped_ref_depth$ref_depth_10x <- rowSums(na.rm=TRUE,white_cases_with_MAF_mutations_reshaped_ref_depth[, 7:ncol(white_cases_with_MAF_mutations_reshaped_ref_depth)] >=10) / nrow(subset(combined_data, (race=='White')))
          white_cases_with_MAF_mutations_reshaped_ref_depth$ref_depth_20x <- rowSums(na.rm=TRUE,white_cases_with_MAF_mutations_reshaped_ref_depth[, 7:ncol(white_cases_with_MAF_mutations_reshaped_ref_depth)] >=20) / nrow(subset(combined_data, (race=='White')))
          white_cases_with_MAF_mutations_reshaped_ref_depth$ref_depth_50x <- rowSums(na.rm=TRUE,white_cases_with_MAF_mutations_reshaped_ref_depth[, 7:ncol(white_cases_with_MAF_mutations_reshaped_ref_depth)] >=50) / nrow(subset(combined_data, (race=='White')))
          white_cases_with_MAF_mutations_reshaped_ref_depth$ref_depth_100x <- rowSums(na.rm=TRUE,white_cases_with_MAF_mutations_reshaped_ref_depth[, 7:ncol(white_cases_with_MAF_mutations_reshaped_ref_depth)] >=100) / nrow(subset(combined_data, (race=='White')))
          
          white_cases_with_MAF_mutations_reshaped_ref_depth <- white_cases_with_MAF_mutations_reshaped_ref_depth[,c(1:6,(ncol(white_cases_with_MAF_mutations_reshaped_ref_depth)-4):ncol(white_cases_with_MAF_mutations_reshaped_ref_depth))]
          
          
          #BASED ON ALT DEPTH
          white_cases_with_MAF_mutations_reshaped_alt_depth <- white_cases_with_MAF_mutations[,c(1,2,3,5,6,23,24,25,28,32)]
          white_cases_with_MAF_mutations_reshaped_alt_depth <- dcast(white_cases_with_MAF_mutations_reshaped_alt_depth, Cancer+race+Chromosome+Start_Position+Tumor_Seq_Ref+Tumor_Seq_Alt ~ submitter_id, value.var = 'tumor_alt_base_count_bam')
          white_cases_with_MAF_mutations_reshaped_alt_depth$alt_depth_1x <- rowSums(na.rm=TRUE,white_cases_with_MAF_mutations_reshaped_alt_depth[, 7:ncol(white_cases_with_MAF_mutations_reshaped_alt_depth)] >=1) / nrow(subset(combined_data, (race=='White')))
          white_cases_with_MAF_mutations_reshaped_alt_depth$alt_depth_5x <- rowSums(na.rm=TRUE,white_cases_with_MAF_mutations_reshaped_alt_depth[, 7:ncol(white_cases_with_MAF_mutations_reshaped_alt_depth)] >=5) / nrow(subset(combined_data, (race=='White')))
          white_cases_with_MAF_mutations_reshaped_alt_depth$alt_depth_10x <- rowSums(na.rm=TRUE,white_cases_with_MAF_mutations_reshaped_alt_depth[, 7:ncol(white_cases_with_MAF_mutations_reshaped_alt_depth)] >=10) / nrow(subset(combined_data, (race=='White')))
          white_cases_with_MAF_mutations_reshaped_alt_depth$alt_depth_20x <- rowSums(na.rm=TRUE,white_cases_with_MAF_mutations_reshaped_alt_depth[, 7:ncol(white_cases_with_MAF_mutations_reshaped_alt_depth)] >=20) / nrow(subset(combined_data, (race=='White')))
          white_cases_with_MAF_mutations_reshaped_alt_depth$alt_depth_50x <- rowSums(na.rm=TRUE,white_cases_with_MAF_mutations_reshaped_alt_depth[, 7:ncol(white_cases_with_MAF_mutations_reshaped_alt_depth)] >=50) / nrow(subset(combined_data, (race=='White')))
          white_cases_with_MAF_mutations_reshaped_alt_depth$alt_depth_100x <- rowSums(na.rm=TRUE,white_cases_with_MAF_mutations_reshaped_alt_depth[, 7:ncol(white_cases_with_MAF_mutations_reshaped_alt_depth)] >=100) / nrow(subset(combined_data, (race=='White')))
          
          #MERGE DATA
          white_cases_with_MAF_mutations_reshaped_alt_depth <- white_cases_with_MAF_mutations_reshaped_alt_depth[,c(1:6,(ncol(white_cases_with_MAF_mutations_reshaped_alt_depth)-4):ncol(white_cases_with_MAF_mutations_reshaped_alt_depth))]
          white_cases_with_MAF_mutations_reshaped <- merge(white_cases_with_MAF_mutations_reshaped_total_depth, white_cases_with_MAF_mutations_reshaped_ref_depth, by=c('Cancer','race','Chromosome','Start_Position','Tumor_Seq_Ref','Tumor_Seq_Alt'))
          white_cases_with_MAF_mutations_reshaped <- merge(white_cases_with_MAF_mutations_reshaped, white_cases_with_MAF_mutations_reshaped_alt_depth, by=c('Cancer','race','Chromosome','Start_Position','Tumor_Seq_Ref','Tumor_Seq_Alt'))
          white_cases_with_MAF_mutations_reshaped <- white_cases_with_MAF_mutations_reshaped[order(white_cases_with_MAF_mutations_reshaped$Start_Position),]
          white_cases_with_MAF_mutations_reshaped[,c(7:ncol(white_cases_with_MAF_mutations_reshaped))] <- round(white_cases_with_MAF_mutations_reshaped[,c(7:ncol(white_cases_with_MAF_mutations_reshaped))]*100, digits=2)
          
          #ADD PERCENTAGE SIGN
          white_cases_with_MAF_mutations_reshaped[,c(7:ncol(white_cases_with_MAF_mutations_reshaped))] <- lapply(white_cases_with_MAF_mutations_reshaped[,c(7:ncol(white_cases_with_MAF_mutations_reshaped))], function(x) paste0(x,"%"))
          rm(white_cases_with_MAF_mutations_reshaped_alt_depth, white_cases_with_MAF_mutations_reshaped_ref_depth, white_cases_with_MAF_mutations_reshaped_total_depth)         }
      
        
        #SUBSET THE POSITIONS SHARED BY BOTH BLACK & WHITE
        #identify positions shared b/w Black and White, and combine with prev data
        shared <- black_cases_with_MAF_mutations_reshaped[black_cases_with_MAF_mutations_reshaped$Start_Position %in% white_cases_with_MAF_mutations_reshaped$Start_Position,]$Start_Position
        
        if (length(shared) > 0){
        Black_White_common_cases_with_MAF_mutations <- rbind((black_cases_with_MAF_mutations_reshaped[black_cases_with_MAF_mutations_reshaped$Start_Position %in% shared,]),(white_cases_with_MAF_mutations_reshaped[white_cases_with_MAF_mutations_reshaped$Start_Position %in% shared,]))
       
        #add cancer name and gene name; also, sort by position
        Black_White_common_cases_with_MAF_mutations$Cancer <- CANCER_NAME
        Black_White_common_cases_with_MAF_mutations$Gene <- GENE
        
        Black_White_common_cases_with_MAF_mutations <- Black_White_common_cases_with_MAF_mutations[order(Black_White_common_cases_with_MAF_mutations$Start_Position),]
        
        #append to dataframe
        Black_White_cases_with_MAF_mutations <- rbind(Black_White_cases_with_MAF_mutations, Black_White_common_cases_with_MAF_mutations)
  
        }
        
        
        #############################################################
        #REARRANGE DATA AND CALCULATE NUMBER/PERCENT OF PATIENTS without DIFFERING COVERAGES: CASES without MAF MUTATIONS
        #############################################################
        
        
        ###
        #BLACK CASES
        ###
        
        #BASED ON TOTAL DEPTH (REF + ALT)
        black_cases_without_MAF_mutations_reshaped_total_depth <- black_cases_without_MAF_mutations[-c(3,4)]
        black_cases_without_MAF_mutations_reshaped_total_depth <- dcast(black_cases_without_MAF_mutations_reshaped_total_depth, Cancer+race+Chromosome+Start_Position+Tumor_Seq_Ref+Tumor_Seq_Alt ~ submitter_id, value.var = 'variant_read_depth_bam_MPILEUP')
        
        if (ncol(black_cases_without_MAF_mutations_reshaped_total_depth) > 7){
          black_cases_without_MAF_mutations_reshaped_total_depth$total_depth_1x <- rowSums(na.rm=TRUE,black_cases_without_MAF_mutations_reshaped_total_depth[, 7:ncol(black_cases_without_MAF_mutations_reshaped_total_depth)] >=1) / nrow(subset(combined_data, (race=='Black')))
          black_cases_without_MAF_mutations_reshaped_total_depth$total_depth_5x <- rowSums(na.rm=TRUE,black_cases_without_MAF_mutations_reshaped_total_depth[, 7:ncol(black_cases_without_MAF_mutations_reshaped_total_depth)] >=5) / nrow(subset(combined_data, (race=='Black')))
          black_cases_without_MAF_mutations_reshaped_total_depth$total_depth_10x <- rowSums(na.rm=TRUE,black_cases_without_MAF_mutations_reshaped_total_depth[, 7:ncol(black_cases_without_MAF_mutations_reshaped_total_depth)] >=10) / nrow(subset(combined_data, (race=='Black')))
          black_cases_without_MAF_mutations_reshaped_total_depth$total_depth_20x <- rowSums(na.rm=TRUE,black_cases_without_MAF_mutations_reshaped_total_depth[, 7:ncol(black_cases_without_MAF_mutations_reshaped_total_depth)] >=20) / nrow(subset(combined_data, (race=='Black')))
          black_cases_without_MAF_mutations_reshaped_total_depth$total_depth_50x <- rowSums(na.rm=TRUE,black_cases_without_MAF_mutations_reshaped_total_depth[, 7:ncol(black_cases_without_MAF_mutations_reshaped_total_depth)] >=50) / nrow(subset(combined_data, (race=='Black')))
          black_cases_without_MAF_mutations_reshaped_total_depth$total_depth_100x <- rowSums(na.rm=TRUE,black_cases_without_MAF_mutations_reshaped_total_depth[, 7:ncol(black_cases_without_MAF_mutations_reshaped_total_depth)] >=100) / nrow(subset(combined_data, (race=='Black')))
          
          black_cases_without_MAF_mutations_reshaped_total_depth <- black_cases_without_MAF_mutations_reshaped_total_depth[,c(1:6,(ncol(black_cases_without_MAF_mutations_reshaped_total_depth)-4):ncol(black_cases_without_MAF_mutations_reshaped_total_depth))]
          
          
          #BASED ON REF DEPTH
          black_cases_without_MAF_mutations_reshaped_ref_depth <- black_cases_without_MAF_mutations[-c(3,4)]
          black_cases_without_MAF_mutations_reshaped_ref_depth <- dcast(black_cases_without_MAF_mutations_reshaped_ref_depth, Cancer+race+Chromosome+Start_Position+Tumor_Seq_Ref+Tumor_Seq_Alt ~ submitter_id, value.var = 'ref_base_count_bam')
          black_cases_without_MAF_mutations_reshaped_ref_depth$ref_depth_1x <- rowSums(na.rm=TRUE,black_cases_without_MAF_mutations_reshaped_ref_depth[, 7:ncol(black_cases_without_MAF_mutations_reshaped_ref_depth)] >=1) / nrow(subset(combined_data, (race=='Black')))
          black_cases_without_MAF_mutations_reshaped_ref_depth$ref_depth_5x <- rowSums(na.rm=TRUE,black_cases_without_MAF_mutations_reshaped_ref_depth[, 7:ncol(black_cases_without_MAF_mutations_reshaped_ref_depth)] >=5) / nrow(subset(combined_data, (race=='Black')))
          black_cases_without_MAF_mutations_reshaped_ref_depth$ref_depth_10x <- rowSums(na.rm=TRUE,black_cases_without_MAF_mutations_reshaped_ref_depth[, 7:ncol(black_cases_without_MAF_mutations_reshaped_ref_depth)] >=10) / nrow(subset(combined_data, (race=='Black')))
          black_cases_without_MAF_mutations_reshaped_ref_depth$ref_depth_20x <- rowSums(na.rm=TRUE,black_cases_without_MAF_mutations_reshaped_ref_depth[, 7:ncol(black_cases_without_MAF_mutations_reshaped_ref_depth)] >=20) / nrow(subset(combined_data, (race=='Black')))
          black_cases_without_MAF_mutations_reshaped_ref_depth$ref_depth_50x <- rowSums(na.rm=TRUE,black_cases_without_MAF_mutations_reshaped_ref_depth[, 7:ncol(black_cases_without_MAF_mutations_reshaped_ref_depth)] >=50) / nrow(subset(combined_data, (race=='Black')))
          black_cases_without_MAF_mutations_reshaped_ref_depth$ref_depth_100x <- rowSums(na.rm=TRUE,black_cases_without_MAF_mutations_reshaped_ref_depth[, 7:ncol(black_cases_without_MAF_mutations_reshaped_ref_depth)] >=100) / nrow(subset(combined_data, (race=='Black')))
          
          black_cases_without_MAF_mutations_reshaped_ref_depth <- black_cases_without_MAF_mutations_reshaped_ref_depth[,c(1:6,(ncol(black_cases_without_MAF_mutations_reshaped_ref_depth)-4):ncol(black_cases_without_MAF_mutations_reshaped_ref_depth))]
          
          
          #BASED ON ALT DEPTH
          black_cases_without_MAF_mutations_reshaped_alt_depth <- black_cases_without_MAF_mutations[-c(3,4)]
          black_cases_without_MAF_mutations_reshaped_alt_depth <- dcast(black_cases_without_MAF_mutations_reshaped_alt_depth, Cancer+race+Chromosome+Start_Position+Tumor_Seq_Ref+Tumor_Seq_Alt ~ submitter_id, value.var = 'alt_base_count_bam')
          black_cases_without_MAF_mutations_reshaped_alt_depth$alt_depth_1x <- rowSums(na.rm=TRUE,black_cases_without_MAF_mutations_reshaped_alt_depth[, 7:ncol(black_cases_without_MAF_mutations_reshaped_alt_depth)] >=1) / nrow(subset(combined_data, (race=='Black')))
          black_cases_without_MAF_mutations_reshaped_alt_depth$alt_depth_5x <- rowSums(na.rm=TRUE,black_cases_without_MAF_mutations_reshaped_alt_depth[, 7:ncol(black_cases_without_MAF_mutations_reshaped_alt_depth)] >=5) / nrow(subset(combined_data, (race=='Black')))
          black_cases_without_MAF_mutations_reshaped_alt_depth$alt_depth_10x <- rowSums(na.rm=TRUE,black_cases_without_MAF_mutations_reshaped_alt_depth[, 7:ncol(black_cases_without_MAF_mutations_reshaped_alt_depth)] >=10) / nrow(subset(combined_data, (race=='Black')))
          black_cases_without_MAF_mutations_reshaped_alt_depth$alt_depth_20x <- rowSums(na.rm=TRUE,black_cases_without_MAF_mutations_reshaped_alt_depth[, 7:ncol(black_cases_without_MAF_mutations_reshaped_alt_depth)] >=20) / nrow(subset(combined_data, (race=='Black')))
          black_cases_without_MAF_mutations_reshaped_alt_depth$alt_depth_50x <- rowSums(na.rm=TRUE,black_cases_without_MAF_mutations_reshaped_alt_depth[, 7:ncol(black_cases_without_MAF_mutations_reshaped_alt_depth)] >=50) / nrow(subset(combined_data, (race=='Black')))
          black_cases_without_MAF_mutations_reshaped_alt_depth$alt_depth_100x <- rowSums(na.rm=TRUE,black_cases_without_MAF_mutations_reshaped_alt_depth[, 7:ncol(black_cases_without_MAF_mutations_reshaped_alt_depth)] >=100) / nrow(subset(combined_data, (race=='Black')))
          
          #MERGE DATA
          black_cases_without_MAF_mutations_reshaped_alt_depth <- black_cases_without_MAF_mutations_reshaped_alt_depth[,c(1:6,(ncol(black_cases_without_MAF_mutations_reshaped_alt_depth)-4):ncol(black_cases_without_MAF_mutations_reshaped_alt_depth))]
          black_cases_without_MAF_mutations_reshaped <- merge(black_cases_without_MAF_mutations_reshaped_total_depth, black_cases_without_MAF_mutations_reshaped_ref_depth, by=c('Cancer','race','Chromosome','Start_Position','Tumor_Seq_Ref','Tumor_Seq_Alt'))
          black_cases_without_MAF_mutations_reshaped <- merge(black_cases_without_MAF_mutations_reshaped, black_cases_without_MAF_mutations_reshaped_alt_depth, by=c('Cancer','race','Chromosome','Start_Position','Tumor_Seq_Ref','Tumor_Seq_Alt'))
          black_cases_without_MAF_mutations_reshaped <- black_cases_without_MAF_mutations_reshaped[order(black_cases_without_MAF_mutations_reshaped$Start_Position),]
          black_cases_without_MAF_mutations_reshaped[,c(7:ncol(black_cases_without_MAF_mutations_reshaped))] <- round(black_cases_without_MAF_mutations_reshaped[,c(7:ncol(black_cases_without_MAF_mutations_reshaped))]*100, digits=2)
          
          #ADD PERCENTAGE SIGN
          black_cases_without_MAF_mutations_reshaped[,c(7:ncol(black_cases_without_MAF_mutations_reshaped))] <- lapply(black_cases_without_MAF_mutations_reshaped[,c(7:ncol(black_cases_without_MAF_mutations_reshaped))], function(x) paste0(x,"%"))
          rm(black_cases_without_MAF_mutations_reshaped_alt_depth, black_cases_without_MAF_mutations_reshaped_ref_depth, black_cases_without_MAF_mutations_reshaped_total_depth)
        }
          ###
          #WHITE CASES
          ###
          
        #BASED ON TOTAL DEPTH (REF + ALT)
        white_cases_without_MAF_mutations_reshaped_total_depth <- white_cases_without_MAF_mutations[-c(3,4)]
        white_cases_without_MAF_mutations_reshaped_total_depth <- dcast(white_cases_without_MAF_mutations_reshaped_total_depth, Cancer+race+Chromosome+Start_Position+Tumor_Seq_Ref+Tumor_Seq_Alt ~ submitter_id, value.var = 'variant_read_depth_bam_MPILEUP')
        
        if (ncol(white_cases_without_MAF_mutations_reshaped_total_depth) > 7){
          white_cases_without_MAF_mutations_reshaped_total_depth$total_depth_1x <- rowSums(na.rm=TRUE,white_cases_without_MAF_mutations_reshaped_total_depth[, 7:ncol(white_cases_without_MAF_mutations_reshaped_total_depth)] >=1) / nrow(subset(combined_data, (race=='White')))
          white_cases_without_MAF_mutations_reshaped_total_depth$total_depth_5x <- rowSums(na.rm=TRUE,white_cases_without_MAF_mutations_reshaped_total_depth[, 7:ncol(white_cases_without_MAF_mutations_reshaped_total_depth)] >=5) / nrow(subset(combined_data, (race=='White')))
          white_cases_without_MAF_mutations_reshaped_total_depth$total_depth_10x <- rowSums(na.rm=TRUE,white_cases_without_MAF_mutations_reshaped_total_depth[, 7:ncol(white_cases_without_MAF_mutations_reshaped_total_depth)] >=10) / nrow(subset(combined_data, (race=='White')))
          white_cases_without_MAF_mutations_reshaped_total_depth$total_depth_20x <- rowSums(na.rm=TRUE,white_cases_without_MAF_mutations_reshaped_total_depth[, 7:ncol(white_cases_without_MAF_mutations_reshaped_total_depth)] >=20) / nrow(subset(combined_data, (race=='White')))
          white_cases_without_MAF_mutations_reshaped_total_depth$total_depth_50x <- rowSums(na.rm=TRUE,white_cases_without_MAF_mutations_reshaped_total_depth[, 7:ncol(white_cases_without_MAF_mutations_reshaped_total_depth)] >=50) / nrow(subset(combined_data, (race=='White')))
          white_cases_without_MAF_mutations_reshaped_total_depth$total_depth_100x <- rowSums(na.rm=TRUE,white_cases_without_MAF_mutations_reshaped_total_depth[, 7:ncol(white_cases_without_MAF_mutations_reshaped_total_depth)] >=100) / nrow(subset(combined_data, (race=='White')))
          
          white_cases_without_MAF_mutations_reshaped_total_depth <- white_cases_without_MAF_mutations_reshaped_total_depth[,c(1:6,(ncol(white_cases_without_MAF_mutations_reshaped_total_depth)-4):ncol(white_cases_without_MAF_mutations_reshaped_total_depth))]
          
          
          #BASED ON REF DEPTH
          white_cases_without_MAF_mutations_reshaped_ref_depth <- white_cases_without_MAF_mutations[-c(3,4)]
          white_cases_without_MAF_mutations_reshaped_ref_depth <- dcast(white_cases_without_MAF_mutations_reshaped_ref_depth, Cancer+race+Chromosome+Start_Position+Tumor_Seq_Ref+Tumor_Seq_Alt ~ submitter_id, value.var = 'ref_base_count_bam')
          white_cases_without_MAF_mutations_reshaped_ref_depth$ref_depth_1x <- rowSums(na.rm=TRUE,white_cases_without_MAF_mutations_reshaped_ref_depth[, 7:ncol(white_cases_without_MAF_mutations_reshaped_ref_depth)] >=1) / nrow(subset(combined_data, (race=='White')))
          white_cases_without_MAF_mutations_reshaped_ref_depth$ref_depth_5x <- rowSums(na.rm=TRUE,white_cases_without_MAF_mutations_reshaped_ref_depth[, 7:ncol(white_cases_without_MAF_mutations_reshaped_ref_depth)] >=5) / nrow(subset(combined_data, (race=='White')))
          white_cases_without_MAF_mutations_reshaped_ref_depth$ref_depth_10x <- rowSums(na.rm=TRUE,white_cases_without_MAF_mutations_reshaped_ref_depth[, 7:ncol(white_cases_without_MAF_mutations_reshaped_ref_depth)] >=10) / nrow(subset(combined_data, (race=='White')))
          white_cases_without_MAF_mutations_reshaped_ref_depth$ref_depth_20x <- rowSums(na.rm=TRUE,white_cases_without_MAF_mutations_reshaped_ref_depth[, 7:ncol(white_cases_without_MAF_mutations_reshaped_ref_depth)] >=20) / nrow(subset(combined_data, (race=='White')))
          white_cases_without_MAF_mutations_reshaped_ref_depth$ref_depth_50x <- rowSums(na.rm=TRUE,white_cases_without_MAF_mutations_reshaped_ref_depth[, 7:ncol(white_cases_without_MAF_mutations_reshaped_ref_depth)] >=50) / nrow(subset(combined_data, (race=='White')))
          white_cases_without_MAF_mutations_reshaped_ref_depth$ref_depth_100x <- rowSums(na.rm=TRUE,white_cases_without_MAF_mutations_reshaped_ref_depth[, 7:ncol(white_cases_without_MAF_mutations_reshaped_ref_depth)] >=100) / nrow(subset(combined_data, (race=='White')))
          
          white_cases_without_MAF_mutations_reshaped_ref_depth <- white_cases_without_MAF_mutations_reshaped_ref_depth[,c(1:6,(ncol(white_cases_without_MAF_mutations_reshaped_ref_depth)-4):ncol(white_cases_without_MAF_mutations_reshaped_ref_depth))]
          
          
          #BASED ON ALT DEPTH
          white_cases_without_MAF_mutations_reshaped_alt_depth <- white_cases_without_MAF_mutations[-c(3,4)]
          white_cases_without_MAF_mutations_reshaped_alt_depth <- dcast(white_cases_without_MAF_mutations_reshaped_alt_depth, Cancer+race+Chromosome+Start_Position+Tumor_Seq_Ref+Tumor_Seq_Alt ~ submitter_id, value.var = 'alt_base_count_bam')
          white_cases_without_MAF_mutations_reshaped_alt_depth$alt_depth_1x <- rowSums(na.rm=TRUE,white_cases_without_MAF_mutations_reshaped_alt_depth[, 7:ncol(white_cases_without_MAF_mutations_reshaped_alt_depth)] >=1) / nrow(subset(combined_data, (race=='White')))
          white_cases_without_MAF_mutations_reshaped_alt_depth$alt_depth_5x <- rowSums(na.rm=TRUE,white_cases_without_MAF_mutations_reshaped_alt_depth[, 7:ncol(white_cases_without_MAF_mutations_reshaped_alt_depth)] >=5) / nrow(subset(combined_data, (race=='White')))
          white_cases_without_MAF_mutations_reshaped_alt_depth$alt_depth_10x <- rowSums(na.rm=TRUE,white_cases_without_MAF_mutations_reshaped_alt_depth[, 7:ncol(white_cases_without_MAF_mutations_reshaped_alt_depth)] >=10) / nrow(subset(combined_data, (race=='White')))
          white_cases_without_MAF_mutations_reshaped_alt_depth$alt_depth_20x <- rowSums(na.rm=TRUE,white_cases_without_MAF_mutations_reshaped_alt_depth[, 7:ncol(white_cases_without_MAF_mutations_reshaped_alt_depth)] >=20) / nrow(subset(combined_data, (race=='White')))
          white_cases_without_MAF_mutations_reshaped_alt_depth$alt_depth_50x <- rowSums(na.rm=TRUE,white_cases_without_MAF_mutations_reshaped_alt_depth[, 7:ncol(white_cases_without_MAF_mutations_reshaped_alt_depth)] >=50) / nrow(subset(combined_data, (race=='White')))
          white_cases_without_MAF_mutations_reshaped_alt_depth$alt_depth_100x <- rowSums(na.rm=TRUE,white_cases_without_MAF_mutations_reshaped_alt_depth[, 7:ncol(white_cases_without_MAF_mutations_reshaped_alt_depth)] >=100) / nrow(subset(combined_data, (race=='White')))
          
          #MERGE DATA
          white_cases_without_MAF_mutations_reshaped_alt_depth <- white_cases_without_MAF_mutations_reshaped_alt_depth[,c(1:6,(ncol(white_cases_without_MAF_mutations_reshaped_alt_depth)-4):ncol(white_cases_without_MAF_mutations_reshaped_alt_depth))]
          white_cases_without_MAF_mutations_reshaped <- merge(white_cases_without_MAF_mutations_reshaped_total_depth, white_cases_without_MAF_mutations_reshaped_ref_depth, by=c('Cancer','race','Chromosome','Start_Position','Tumor_Seq_Ref','Tumor_Seq_Alt'))
          white_cases_without_MAF_mutations_reshaped <- merge(white_cases_without_MAF_mutations_reshaped, white_cases_without_MAF_mutations_reshaped_alt_depth, by=c('Cancer','race','Chromosome','Start_Position','Tumor_Seq_Ref','Tumor_Seq_Alt'))
          white_cases_without_MAF_mutations_reshaped <- white_cases_without_MAF_mutations_reshaped[order(white_cases_without_MAF_mutations_reshaped$Start_Position),]
          white_cases_without_MAF_mutations_reshaped[,c(7:ncol(white_cases_without_MAF_mutations_reshaped))] <- round(white_cases_without_MAF_mutations_reshaped[,c(7:ncol(white_cases_without_MAF_mutations_reshaped))]*100, digits=2)
          
          #ADD PERCENTAGE SIGN
          white_cases_without_MAF_mutations_reshaped[,c(7:ncol(white_cases_without_MAF_mutations_reshaped))] <- lapply(white_cases_without_MAF_mutations_reshaped[,c(7:ncol(white_cases_without_MAF_mutations_reshaped))], function(x) paste0(x,"%"))
          rm(white_cases_without_MAF_mutations_reshaped_alt_depth, white_cases_without_MAF_mutations_reshaped_ref_depth, white_cases_without_MAF_mutations_reshaped_total_depth)
        }
         
          
          
          #SUBSET THE POSITIONS SHARED BY BOTH BLACK & WHITE
          #identify positions shared b/w Black and White, and combine without prev data
          shared <- black_cases_without_MAF_mutations_reshaped[black_cases_without_MAF_mutations_reshaped$Start_Position %in% white_cases_without_MAF_mutations_reshaped$Start_Position,]$Start_Position
          
          if (length(shared) > 0){
            Black_White_common_cases_without_MAF_mutations <- rbind((black_cases_without_MAF_mutations_reshaped[black_cases_without_MAF_mutations_reshaped$Start_Position %in% shared,]),(white_cases_without_MAF_mutations_reshaped[white_cases_without_MAF_mutations_reshaped$Start_Position %in% shared,]))
            
            #add cancer name and gene name; also, sort by position
            Black_White_common_cases_without_MAF_mutations$Cancer <- CANCER_NAME
            Black_White_common_cases_without_MAF_mutations$Gene <- GENE
            
            Black_White_common_cases_without_MAF_mutations <- Black_White_common_cases_without_MAF_mutations[order(Black_White_common_cases_without_MAF_mutations$Start_Position),]
            #append to dataframe
            Black_White_cases_without_MAF_mutations <- rbind(Black_White_cases_without_MAF_mutations, Black_White_common_cases_without_MAF_mutations)
          }
      }
  }      
      
      
      
    
#REARRANGE
Black_White_cases_with_MAF_mutations <- Black_White_cases_with_MAF_mutations[,c(22,1:21)]
Black_White_cases_without_MAF_mutations <- Black_White_cases_without_MAF_mutations[,c(22,1:21)]

rm(list=setdiff(ls(), c("Black_White_cases_with_MAF_mutations","Black_White_cases_without_MAF_mutations","GENE")))
write.table(x=Black_White_cases_with_MAF_mutations, file=paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/racial_comparison_loci_tables/with_MAF_mutations_', GENE),row.names = FALSE, sep='\t')
write.table(x=Black_White_cases_without_MAF_mutations, file=paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/racial_comparison_loci_tables/without_MAF_mutations_', GENE),row.names = FALSE, sep='\t')
