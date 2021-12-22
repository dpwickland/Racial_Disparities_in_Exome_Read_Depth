#install packages using jupyter conda environment
library(forcats)
library(ggplot2)
library(scales)
library(reshape2)
library(graphics)
library(ggpubr)
library(tidyr)
library(stringr)
library(plyr)
library(lemon)
library(dplyr)
library(ggpubr)
library(shades)
library(gridExtra)
library(grid)
library(gridGraphics)
library(cowplot)


theme_set(theme_grey())


######################################
#MASTER TABLE
######################################

PLOT_BY_RACE_AND_SAMPLING_SITE <- function(MAPPED_READS_OR_ALL_READS,CANCER_TYPES,EXOME_OR_RNASEQ){
  
  if (CANCER_TYPES=='all'){
    CANCER_LIST <- c('BRCA','UCEC','PRAD','LUAD','COAD','KIRP','KIRC')
    
    #Create dataframe to hold master table for all cancers
    master_table_master <- data.frame()
    
    if(MAPPED_READS_OR_ALL_READS != 'mapped_BRCA_reads_within_nimblegen_capture'){
      
      #Load data
      for (CANCER_NAME in CANCER_LIST){
        
        #load overall depths data
        master_table <- read.table(paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/overall_read_depths_and_sample_information/',CANCER_NAME,'_',MAPPED_READS_OR_ALL_READS,'_',EXOME_OR_RNASEQ,'_master_table_final.txt'),header=TRUE,sep='\t')
        master_table_master <- rbind(master_table_master,master_table)
      }
    }
  }
  
  if (CANCER_TYPES != 'all'){
    CANCER_NAME <- CANCER_TYPES
    master_table_master <- read.table(paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/overall_read_depths_and_sample_information/',CANCER_NAME,'_',MAPPED_READS_OR_ALL_READS,'_',EXOME_OR_RNASEQ,'_master_table_final.txt'),header=TRUE,sep='\t')
    
    
  }
  
  
  #remove entries that have no exome data; some in dataframes have RNA-seq but no exome
  master_table_master <- subset(master_table_master, (bam_full_path != 'NA'))
  
  #set order for race 
  master_table_master <- subset(master_table_master, (race_PCA == 'White' | race_PCA == 'Black'))
  master_table_master$race_PCA <- factor(master_table_master$race_PCA, levels = c("White","Black"),ordered=TRUE)
  
  #rename tissue type and set order
  master_table_master$tissue.definition <- sub("Blood_Normal","Patient Germline",master_table_master$tissue.definition)
  master_table_master$tissue.definition <- sub("Primary_Tumor","Primary Tumor",master_table_master$tissue.definition)
  master_table_master$tissue.definition <- factor(master_table_master$tissue.definition, levels=c('Primary Tumor','Patient Germline'))
  master_table_master$Cancer <- factor(master_table_master$Cancer,levels=c('BRCA','UCEC','PRAD','LUAD','COAD','KIRP','KIRC'),ordered=TRUE)
  
  #if any capture kit listed multiple times for a patient, combine them into single row, so that no patient is counted more than once for same tissue type
  CAPTURE_KIT_FIXED <- master_table_master %>% group_by(submitter_id_long,tissue.definition) %>% dplyr::summarise(CAPTURE_KIT_FIXED = toString(CAPTURE_KIT))
  CAPTURE_KIT_FIXED$CAPTURE_KIT_FIXED <- gsub(", ","|",CAPTURE_KIT_FIXED$CAPTURE_KIT_FIXED)
  
  master_table_master <- merge(master_table_master, CAPTURE_KIT_FIXED, by=c('submitter_id_long','tissue.definition'))
  master_table_master <- master_table_master[(names(master_table_master)!='CAPTURE_KIT')]
  
  master_table_master <- unique(master_table_master)
  
  #simplify capture kit
    master_table_master$CAPTURE_KIT_simp <-  ifelse(grepl('hg18',master_table_master$CAPTURE_KIT_FIXED),'NimbleGen hg18 Exome v2',ifelse(grepl('v3.0',master_table_master$CAPTURE_KIT_FIXED),'NimbleGen SeqCap EZ Exome v3',ifelse(grepl('v2.0',master_table_master$CAPTURE_KIT_FIXED),'NimbleGen SeqCap EZ Exome v2',ifelse(grepl('V2.0',master_table_master$CAPTURE_KIT_FIXED),'NimbleGen SeqCap EZ Exome v2',ifelse(grepl('Custom',master_table_master$CAPTURE_KIT_FIXED),'Custom V2 Exome Bait',ifelse(grepl('Gapfiller',master_table_master$CAPTURE_KIT_FIXED),'Gapfiller_7m',ifelse(grepl('Rome',master_table_master$CAPTURE_KIT_FIXED),'Roche SeqCap EZ HGSC VCRome',ifelse(grepl('SureSelect',master_table_master$CAPTURE_KIT_FIXED),'Agilent SureSelect Human All Exon v2','Not reported')))))))) 
  
  
  #clean up capture kit and order according to date
  master_table_master$CAPTURE_KIT_simp <- gsub("SureSelect Human All Exon 38 Mb v2","Agilent SureSelect Human All Exon v2",master_table_master$CAPTURE_KIT_simp)
  
  
  #reorder capture kit simp
  master_table_master$CAPTURE_KIT_simp <- factor(master_table_master$CAPTURE_KIT_simp, levels=c('Agilent SureSelect Human All Exon v2','Custom V2 Exome Bait','NimbleGen hg18 Exome v2','NimbleGen SeqCap EZ Exome v2','NimbleGen SeqCap EZ Exome v3','Roche SeqCap EZ HGSC VCRome','Gapfiller_7m','Not reported'),ordered=TRUE)
  
  return(master_table_master)
}


ALL_READS_FOR_PLOTTING <- unique(PLOT_BY_RACE_AND_SAMPLING_SITE('all_reads','all','exome'))
MAPPED_READS_FOR_PLOTTING <- unique(PLOT_BY_RACE_AND_SAMPLING_SITE('mapped_reads','all','exome'))
UNMAPPED_READS_FOR_PLOTTING <- unique(PLOT_BY_RACE_AND_SAMPLING_SITE('unmapped_reads','all','exome'))

SeqCapV2_BRCA_READS_FOR_PLOTTING <- unique(PLOT_BY_RACE_AND_SAMPLING_SITE('mapped_reads_within_nimblegen_SeqCap_EZ_Exome_v2_hg38liftover_SORTED','BRCA','exome'))
SeqCapV3_BRCA_READS_FOR_PLOTTING <- unique(PLOT_BY_RACE_AND_SAMPLING_SITE('mapped_reads_within_nimblegen_SeqCap_EZ_Exome_v3_hg19_capture_targets_hg38liftover_SORTED','BRCA','exome'))
hg18_BRCA_READS_FOR_PLOTTING <- unique(PLOT_BY_RACE_AND_SAMPLING_SITE('mapped_reads_within_nimblegen_hg18_nimblegen_exome_version_2_hg38liftover_SORTED','BRCA','exome'))


#############################
#ADD SEQUENCING DATE TO TABLE
############################
#https://api.gdc.cancer.gov/files?filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22!=%22%2C%22content%22%3A%7B%22field%22%3A%22cases.submitter_id%22%2C%22value%22%3A%5B%22TCGA-4V-MIUM%22%5D%7D%7D%2C%7B%22op%22%3A%22%3D%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_type%22%2C%22value%22%3A%22Aligned%20Reads%22%7D%7D%5D%7D&format=tsv&fields=file_id,file_name,cases.submitter_id,cases.case_id,data_category,data_type,cases.samples.tumor_descriptor,cases.samples.tissue_type,cases.samples.sample_type,cases.samples.submitter_id,cases.samples.sample_id,analysis.workflow_type,cases.project.project_id,cases.samples.portions.analytes.aliquots.aliquot_id,cases.samples.portions.analytes.aliquots.submitter_id,analysis.metadata.read_groups.target_capture_kit_name,analysis.metadata.read_groups.target_capture_kit_target_region,analysis.metadata.read_groups.target_capture_kit_name&size=100000

SEQ_DATE <- read.table('/home/mayo/m187735/s212975.Wickland_Immunomics/TCGA_metadata/seq_date.txt',sep='\t',header=TRUE)
SEQ_DATE <- subset(SEQ_DATE, Workflow_Type=='BWA with Mark Duplicates and Cocleaning')[,c('submitter_id_long','Sequencing_Date')]
SEQ_DATE$Sequencing_Date <- sub("^([^-]*).*", "\\1",SEQ_DATE$Sequencing_Date)

MAPPED_READS_FOR_PLOTTING <- merge(MAPPED_READS_FOR_PLOTTING, SEQ_DATE, by='submitter_id_long',all.x=TRUE)


##################################
#BRCA CAPTURE KITS -- KEEP ONLY THE SAMPLES CAPTURED BY THEIR RESPECTIVE KITS
##################################

#if sample called by two kits, just consider as having been captured by smallest
table(SeqCapV2_BRCA_READS_FOR_PLOTTING$CAPTURE_KIT_simp)
SeqCapV2_BRCA_READS_FOR_PLOTTING <- subset(SeqCapV2_BRCA_READS_FOR_PLOTTING, CAPTURE_KIT_simp == 'NimbleGen SeqCap EZ Exome v2')
SeqCapV3_BRCA_READS_FOR_PLOTTING <- subset(SeqCapV3_BRCA_READS_FOR_PLOTTING, CAPTURE_KIT_simp == 'NimbleGen SeqCap EZ Exome v3')
hg18_BRCA_READS_FOR_PLOTTING <- subset(hg18_BRCA_READS_FOR_PLOTTING, CAPTURE_KIT_simp == 'NimbleGen hg18 Exome v2')

BRCA_three_nimblegen_kits_overall_reads <- rbind(hg18_BRCA_READS_FOR_PLOTTING, SeqCapV2_BRCA_READS_FOR_PLOTTING,SeqCapV3_BRCA_READS_FOR_PLOTTING)

BRCA_three_nimblegen_kits_overall_reads <- merge(BRCA_three_nimblegen_kits_overall_reads, SEQ_DATE, by='submitter_id_long',all.x=TRUE)


##################################
#CAPTURE KIT TABLE -
##################################
#if sample called by two kits, just consider as having been captured by smallest
CAPTURE_TABLE <- (droplevels(MAPPED_READS_FOR_PLOTTING[,c('race_PCA','submitter_id','tissue.definition','Cancer','CAPTURE_KIT_FIXED')]) %>%
                    dplyr::group_by(Cancer,tissue.definition,race_PCA,CAPTURE_KIT_FIXED,.drop=FALSE) %>%
                    dplyr::summarise(number_of_samples=dplyr::n()) ) 
CAPTURE_TABLE <- subset(CAPTURE_TABLE, number_of_samples != 0)

CAPTURE_TABLE <- CAPTURE_TABLE %>% spread(tissue.definition,number_of_samples)


CAPTURE_TABLE_WHITE <- subset(CAPTURE_TABLE, race_PCA=='White')
names(CAPTURE_TABLE_WHITE)[4:5] <-c( "Primary Tumor White","Patient Germline White")

CAPTURE_TABLE_BLACK <- subset(CAPTURE_TABLE, race_PCA=='Black')
names(CAPTURE_TABLE_BLACK)[4:5] <-c( "Primary Tumor Black","Patient Germline Black")

CAPTURE_TABLE <- merge(CAPTURE_TABLE_WHITE[-c(2)],CAPTURE_TABLE_BLACK[-c(2)],by=c('Cancer','CAPTURE_KIT_FIXED'),all=TRUE)

CAPTURE_TABLE$CAPTURE_KIT_simp <-  ifelse(grepl('hg18',CAPTURE_TABLE$CAPTURE_KIT_FIXED),'NimbleGen hg18 Exome v2',ifelse(grepl('v3.0',CAPTURE_TABLE$CAPTURE_KIT_FIXED),'NimbleGen SeqCap EZ Exome v3',ifelse(grepl('v2.0',CAPTURE_TABLE$CAPTURE_KIT_FIXED),'NimbleGen SeqCap EZ Exome v2',ifelse(grepl('V2.0',CAPTURE_TABLE$CAPTURE_KIT_FIXED),'NimbleGen SeqCap EZ Exome v2',ifelse(grepl('Custom',CAPTURE_TABLE$CAPTURE_KIT_FIXED),'Custom V2 Exome Bait',ifelse(grepl('Gapfiller',CAPTURE_TABLE$CAPTURE_KIT_FIXED),'Gapfiller_7m',ifelse(grepl('Rome',CAPTURE_TABLE$CAPTURE_KIT_FIXED),'Roche SeqCap EZ HGSC VCRome',ifelse(grepl('SureSelect',CAPTURE_TABLE$CAPTURE_KIT_FIXED),'Agilent SureSelect Human All Exon v2','Not reported'))))))))  

CAPTURE_TABLE$CAPTURE_KIT_simp <- factor(CAPTURE_TABLE$CAPTURE_KIT_simp, levels=c('Agilent SureSelect Human All Exon v2','Custom V2 Exome Bait','NimbleGen hg18 Exome v2','NimbleGen SeqCap EZ Exome v2','NimbleGen SeqCap EZ Exome v3','Roche SeqCap EZ HGSC VCRome','Not reported'),ordered=TRUE)






######################################
#READ COUNT PER PATIENT (BAM)
######################################
PLOT_OVERALL_READ_DEPTHS_BY_RACE <- function(MAPPED_READS_OR_ALL_READS,FIG_NUMBER){
  
  if (MAPPED_READS_OR_ALL_READS=='all_reads'){DF_TO_PLOT <- ALL_READS_FOR_PLOTTING; TITLE <- 'Total number of reads per exome'; UNITS <- 'Millions'; N_POSITION=-50000000}
  if (MAPPED_READS_OR_ALL_READS=='mapped_reads'){DF_TO_PLOT <- MAPPED_READS_FOR_PLOTTING; TITLE <- 'Number of mapped reads per exome'; UNITS <- 'Millions'; N_POSITION=-50000000} 
  if (MAPPED_READS_OR_ALL_READS=='unmapped_reads'){DF_TO_PLOT <- UNMAPPED_READS_FOR_PLOTTING; TITLE <- 'Number of unmapped reads per exome'; UNITS <- 'Thousands'; N_POSITION=10} 
  # if (MAPPED_READS_OR_ALL_READS=='mapped_HLA_reads'){DF_TO_PLOT <- MAPPED_HLA_READS_FOR_PLOTTING; TITLE <- '# of mapped HLA-I and HLA-II reads per patient'; UNITS <- 'thousands'; N_POSITION=-200000} 
  
  
  DF_TO_PLOT <- subset(DF_TO_PLOT, (read_depth > 0))
  
  #calculate summary statistics
  Counts <- aggregate(read_depth ~ tissue.definition + Cancer + race_PCA, DF_TO_PLOT, length)
  Means <- aggregate(read_depth ~ tissue.definition + Cancer + race_PCA, DF_TO_PLOT, mean)
  if(MAPPED_READS_OR_ALL_READS=='all_reads' | MAPPED_READS_OR_ALL_READS=='mapped_reads'){Means$read_depth <- sprintf("%1.1f",round(Means$read_depth,1)/1000000); COORD_TOP=7200000;COORD_BOT=3990000}
  if(MAPPED_READS_OR_ALL_READS=='unmapped_reads'){Means$read_depth <- sprintf("%1.1f",round(Means$read_depth,1)/1000);COORD_TOP=300;COORD_BOT=61}
  
  
  #create labels for the annotated text below each facet row
  tumor_label_counts <- data.frame(race_PCA='White',tissue.definition=c('Primary Tumor'), read_depth=c(14000000),Cancer=c('BRCA','UCEC','PRAD','LUAD','COAD','KIRP','KIRC'),label=c('n','','','','','',''))
  tumor_label_counts$Cancer <- factor(tumor_label_counts$Cancer, levels=c('BRCA','UCEC','PRAD','LUAD','COAD','KIRP','KIRC'))
  germline_label_counts <- data.frame(race_PCA='White',tissue.definition=c('Patient Germline'), read_depth=c(14000000),Cancer=c('BRCA','UCEC','PRAD','LUAD','COAD','KIRP','KIRC'),label=c('n','','','','','',''))
  germline_label_counts$Cancer <- factor(germline_label_counts$Cancer, levels=c('BRCA','UCEC','PRAD','LUAD','COAD','KIRP','KIRC'))
  
  tumor_label_mean <- data.frame(race_PCA='White',tissue.definition=c('Primary Tumor'), read_depth=c(14000000),Cancer=c('BRCA','UCEC','PRAD','LUAD','COAD','KIRP','KIRC'),label=c('bar(x)','','','','','',''))
  tumor_label_mean$Cancer <- factor(tumor_label_mean$Cancer, levels=c('BRCA','UCEC','PRAD','LUAD','COAD','KIRP','KIRC'))
  germline_label_mean <- data.frame(race_PCA='White',tissue.definition=c('Patient Germline'), read_depth=c(14000000),Cancer=c('BRCA','UCEC','PRAD','LUAD','COAD','KIRP','KIRC'),label=c('bar(x)','','','','','',''))
  germline_label_mean$Cancer <- factor(germline_label_mean$Cancer, levels=c('BRCA','UCEC','PRAD','LUAD','COAD','KIRP','KIRC'))
  
  DF_TO_PLOT$Cancer <- factor(DF_TO_PLOT$Cancer, levels=c('BRCA','UCEC','PRAD','LUAD','COAD','KIRP','KIRC'))
  
  all_reads_by_race_PCA <- ggplot(DF_TO_PLOT, aes(x=factor(race_PCA),y=read_depth,fill=race_PCA)) + xlab("") + ylab(paste0(UNITS, " of reads per exome (log2 scale)")) + facet_grid(tissue.definition~Cancer,labeller = label_wrap_gen(width=21)) + labs(title=TITLE) + 
    geom_violin() + geom_sina(size=0.25) +  stat_summary(geom='crossbar',fun='median',color='white',fatten=2,width=0.35) + guides(fill = guide_legend(override.aes = list(linetype = 0))) + 
    theme(strip.text.x = element_text(size = 24),#,hjust=0),
          strip.text.y = element_text(size = 24),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size=22,margin=margin(r=-1)),
          axis.title.y = element_text(size=24,margin=margin(t=0,r=10,b=0,l=0)),
          axis.ticks.x = element_blank(), 
          panel.spacing.y=unit(5,"lines"), #to control distance betweeen facet rows
          #     plot.tag = element_text(size = 30,face='bold'),
          #    plot.tag.position = c(-0.01, 0.98),
          strip.background = element_blank(),
          panel.background = element_rect(color = 'black',fill=NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust=0.5,size=30,face='bold',margin=margin(t=0,r=0,b=15,l=0)),
          legend.position = 'bottom',legend.direction='horizontal',
          legend.text=element_text(size=34),
          legend.background = element_rect(linetype='solid', colour='black'),
          legend.title=element_blank()) + 
    #geom_point() +
    
    
    #add text annotations
    geom_text(data = subset(Counts, tissue.definition == "Primary Tumor"),aes(label=read_depth, y=COORD_TOP,vjust=0.2),size=7) + 
    geom_text(data = subset(Means, tissue.definition == "Primary Tumor"),aes(label=read_depth, y=COORD_BOT,vjust=0.2),size=7) + 
    geom_text(data= tumor_label_counts,aes(label=label,y=COORD_TOP,vjust=0.2,hjust=4.3,fontface="bold"),size=7) +
    geom_text(data= tumor_label_mean,aes(label=label,y=COORD_BOT,vjust=0.2,hjust=5.15,fontface="bold"),size=7,parse = TRUE) +
    geom_text(data= tumor_label_mean,aes(label=label,y=COORD_BOT,vjust=0.2,hjust=5.15,fontface="bold"),size=7,parse=TRUE) + #}} +
    geom_text(data= tumor_label_mean,aes(label=label,y=COORD_BOT,vjust=0.2,hjust=5.15,fontface="bold"),size=7.05,parse=TRUE) + #}} +
    geom_text(data= tumor_label_mean,aes(label=label,y=COORD_BOT,vjust=0.2,hjust=5.15,fontface="bold"),size=7.10,parse=TRUE) + #}} +
    geom_text(data= tumor_label_mean,aes(label=label,y=COORD_BOT,vjust=0.2,hjust=5.15,fontface="bold"),size=7.12,parse=TRUE) + #}} +
    geom_text(data= tumor_label_mean,aes(label=label,y=COORD_BOT,vjust=0.2,hjust=5.15,fontface="bold"),size=7.15,parse=TRUE) + #}} +
    geom_text(data= tumor_label_mean,aes(label=label,y=COORD_BOT,vjust=0.2,hjust=5.15,fontface="bold"),size=7.20,parse=TRUE) + #}} +
    geom_text(data= tumor_label_mean,aes(label=label,y=COORD_BOT,vjust=0.2,hjust=5.15,fontface="bold"),size=7.18,parse=TRUE) + #}} +
    
    
    #  {if(EXOME_OR_RNASEQ=='exome'){
    geom_text(data = subset(Counts, tissue.definition == "Patient Germline"),aes(label=read_depth, y=COORD_TOP,vjust=0.2),size=7) + 
    geom_text(data = subset(Means, tissue.definition == "Patient Germline"),aes(label=read_depth, y=COORD_BOT,vjust=0.2),size=7) + 
    geom_text(data= germline_label_counts,aes(label=label,y=COORD_TOP,vjust=0.2,hjust=4.3,fontface="bold"),size=7) +
    geom_text(data= germline_label_mean,aes(label=label,y=COORD_BOT,vjust=0.2,hjust=5.15,fontface="bold"),size=7,parse=TRUE) + #}} +
    geom_text(data= germline_label_mean,aes(label=label,y=COORD_BOT,vjust=0.2,hjust=5.15,fontface="bold"),size=7.05,parse=TRUE) + #}} +
    geom_text(data= germline_label_mean,aes(label=label,y=COORD_BOT,vjust=0.2,hjust=5.15,fontface="bold"),size=7.10,parse=TRUE) + #}} +
    geom_text(data= germline_label_mean,aes(label=label,y=COORD_BOT,vjust=0.2,hjust=5.15,fontface="bold"),size=7.12,parse=TRUE) + #}} +
    geom_text(data= germline_label_mean,aes(label=label,y=COORD_BOT,vjust=0.2,hjust=5.15,fontface="bold"),size=7.15,parse=TRUE) + #}} +
    geom_text(data= germline_label_mean,aes(label=label,y=COORD_BOT,vjust=0.2,hjust=5.15,fontface="bold"),size=7.20,parse=TRUE) + #}} +
    geom_text(data= germline_label_mean,aes(label=label,y=COORD_BOT,vjust=0.2,hjust=5.15,fontface="bold"),size=7.18,parse=TRUE) + #}} +
    
    
    
    #allow table outside plot margins
    {if(MAPPED_READS_OR_ALL_READS=='all_reads' | MAPPED_READS_OR_ALL_READS=='mapped_reads') { coord_cartesian(ylim=c(12000000,400000000),clip='off')}}    +  
    {if(MAPPED_READS_OR_ALL_READS=='unmapped_reads') {coord_cartesian(ylim=c(5500,35000000),clip='off') }}    +  
    
    
    scale_fill_manual(values=c('#EF5350','#42A5F5')) + #MANUALLY SET FILL COLORS  
    
    
    #means comparisons
    stat_compare_means(size=7,method='wilcox.test',mapping=aes(label=..p.signif..),
                       comparisons = list(c("Black","White")),bracket.size = NA,tip.length = 0,color='black',
                       #label.y=400000000, #adjust position of p-value
                       symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1,100), symbols = c("p<0.0001", "p<0.001", "p<0.01", "p<0.05", "",""))) + 
    
    #   {if(MAPPED_READS_OR_ALL_READS=='all_reads' | MAPPED_READS_OR_ALL_READS=='mapped_reads') {expand_limits(y=-500000)}}    +  
    #    {if(MAPPED_READS_OR_ALL_READS=='all_reads' | MAPPED_READS_OR_ALL_READS=='mapped_reads') {labs(tag = "A.")}}    +  
    #   {if(MAPPED_READS_OR_ALL_READS=='unmapped_reads') {labs(tag = "B.")}}    +  
    
    {if(MAPPED_READS_OR_ALL_READS=='all_reads' | MAPPED_READS_OR_ALL_READS=='mapped_reads') {scale_y_continuous(expand = expansion(mult = c(0, 0.30)),labels=unit_format(unit = "",scale = 1e-6,accuracy=1,big.mark = ','),trans=pseudo_log_trans(base=2),breaks = c(25000000,50000000,100000000,200000000,400000000))}
      
      else if (MAPPED_READS_OR_ALL_READS=='unmapped_reads') {scale_y_continuous(expand = expansion(mult = c(.15, 0.30)),labels=unit_format(unit = "",scale = 1e-3,accuracy=1,big.mark = ','),trans=pseudo_log_trans(base=2),breaks = c(5000,25000,100000,500000,2000000,10000000))}
      
      else if (MAPPED_READS_OR_ALL_READS=='mapped_HLA_reads') {scale_y_continuous(labels = unit_format(unit = "", scale = 1e-3), expand = expansion(mult = c(.1, 0.2)),breaks = 2^seq(1,12, by = 2))}}
  
  
}

OVERALL_READ_DEPTHS_BY_RACE_ALL_READS <- PLOT_OVERALL_READ_DEPTHS_BY_RACE('all_reads','S1A')
OVERALL_READ_DEPTHS_BY_RACE_MAPPED_READS <- PLOT_OVERALL_READ_DEPTHS_BY_RACE('mapped_reads','1A')
OVERALL_READ_DEPTHS_BY_RACE_UNMAPPED_READS <- PLOT_OVERALL_READ_DEPTHS_BY_RACE('unmapped_reads','S1B')
#PLOT_OVERALL_READ_DEPTHS_BY_RACE('mapped_HLA_reads')





#####################################
#TUMOR PURITY
######################################

Counts_TP <- aggregate(Tumor_Purity ~ tissue.definition + Cancer + race_PCA, MAPPED_READS_FOR_PLOTTING, length)

MAPPED_READS_FOR_PLOTTING$Facet <- '' #add empty facet strip to line up plots

Tumor_Purity <- ggplot(MAPPED_READS_FOR_PLOTTING, aes(fill=race_PCA,x=race_PCA, y=as.numeric(as.character(Tumor_Purity)))) +   facet_grid(Facet~Cancer) + theme(strip.text = element_text(face="bold", size=8)) + 
   geom_violin() + geom_sina(size=0.25) +  stat_summary(geom='crossbar',fun='median',color='white',fatten=2,width=0.35) + guides(fill = guide_legend(override.aes = list(linetype = 0))) + 
  #+ geom_jitter(size=0.5,width=0.025) +
  ylab("Tumor purity") +  xlab("") + labs(title = ("Tumor purity per exome")) +  
  #geom_point() + 
  
  theme(strip.text.x = element_text(size = 24,face='plain'),#,hjust=0),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=22,margin=margin(t=0,r=5,b=0,l=0)),
        axis.title.y = element_text(size=24,margin=margin(t=0,r=10,b=0,l=0)),
        axis.ticks.x = element_blank(), 
        plot.title = element_text(hjust=0.5,size=30,face='bold',margin=margin(t=0,r=0,b=15,l=0)),
        #  plot.tag = element_text(size = 30,face='bold'),
        # plot.tag.position = c(-0.01, 0.98),
        strip.background = element_blank(),
        panel.background = element_rect(color = 'black',fill=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(2,0,0,0), "lines"), #top, right, bottom, left
        legend.position = 'bottom',legend.direction='horizontal',
        legend.background = element_rect(linetype='solid', colour='black'),
        legend.text=element_text(size=34),        
        legend.title=element_blank()) + #labs(tag = "C.") +
  
  
  
  #coord_cartesian(ylim = c(0,1.35))+
  scale_fill_manual(values=c('#EF5350','#42A5F5')) + #MANUALLY SET FILL COLORS
  scale_color_manual(values='orange') + 
  #means comparisons
  stat_compare_means(size=7,method='wilcox.test',mapping=aes(label=..p.signif..),
                     comparisons = list(c("Black","White")),bracket.size = NA,tip.length = 0,color='black',
                     #label.y=400000000, #adjust position of p-value
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("p<0.0001", "p<0.001", "p<0.01", "p<0.05", ""))) +
  expand_limits(y=-0.065) + scale_y_continuous(breaks=c(0,0.50,1),expand = expansion(mult = c(0.07, 0.30))) + 
  
  #add text annotation
  geom_text(data = subset(Counts_TP, tissue.definition == "Primary Tumor"),aes(label=Tumor_Purity, y=-0.32,vjust=0),size=7)   +
  coord_cartesian(ylim=c(0,1.2),clip='off')


#OUTPUT FIGURE 1
pdf(paste0("/home/mayo/m187735/s212975.Wickland_Immunomics/GitHub/Racial_Disparities_in_Exome_Read_Depth/MS_plots/FIG1.pdf"),height=14,width=19,onefile = FALSE)
ggarrange(NULL,OVERALL_READ_DEPTHS_BY_RACE_MAPPED_READS,NULL,NULL,NULL,Tumor_Purity,NULL,NULL,nrow=4,ncol=2,common.legend = TRUE,legend = 'bottom',align='hv',widths=c(0.04,1),heights=c(0.9,0.07,0.50,0.11),labels=c('A','','','','B','','',''),font.label = list(size = 32, color = "black"))  + theme(plot.margin = margin(b=0.5,l=0.1,t=0.1,r=0.1, "cm")) #incr marrgin at bottom
dev.off()

#png(paste0("/home/mayo/m187735/s212975.Wickland_Immunomics/GitHub/Racial_Disparities_in_Exome_Read_Depth/MS_plots/FIG1.png"),height=13.5,width=19,res=500,units='in')
#ggarrange(NULL,OVERALL_READ_DEPTHS_BY_RACE_MAPPED_READS,NULL,NULL,NULL,Tumor_Purity,NULL,NULL,nrow=4,ncol=2,common.legend = TRUE,legend = 'bottom',align='hv',widths=c(0.04,1),heights=c(0.9,0.07,0.50,0.11),labels=c('A','','','','B','','',''),font.label = list(size = 32, color = "black"))  + theme(plot.margin = margin(b=0.5,l=0.1,t=0.1,r=0.1, "cm")) #incr marrgin at bottom
#dev.off()

#OUTPUT FIGURE S2
#pdf(paste0("/home/mayo/m187735/s212975.Wickland_Immunomics/GitHub/Racial_Disparities_in_Exome_Read_Depth/MS_plots/FIGS2.pdf"),height=17.5,width=19,onefile = FALSE)
#ggarrange(NULL,OVERALL_READ_DEPTHS_BY_RACE_ALL_READS,NULL,NULL,NULL,OVERALL_READ_DEPTHS_BY_RACE_UNMAPPED_READS,NULL,NULL,nrow=4,ncol=2,common.legend = TRUE,legend = 'bottom',align='hv',widths=c(0.04,1),heights=c(0.9,0.12,0.90,0.11),labels=c('A','','','','B','','',''),font.label = list(size = 32, color = "black"))  + theme(plot.margin = margin(b=0.5,l=0.1,t=0.1,r=0.1, "cm")) #incr marrgin at bottom
#dev.off()

png(paste0("/home/mayo/m187735/s212975.Wickland_Immunomics/GitHub/Racial_Disparities_in_Exome_Read_Depth/MS_plots/FIGS2.png"),height=17.5,width=19,res=500,units='in')
ggarrange(NULL,OVERALL_READ_DEPTHS_BY_RACE_ALL_READS,NULL,NULL,NULL,OVERALL_READ_DEPTHS_BY_RACE_UNMAPPED_READS,NULL,NULL,nrow=4,ncol=2,common.legend = TRUE,legend = 'bottom',align='hv',widths=c(0.04,1),heights=c(0.9,0.12,0.90,0.11),labels=c('A','','','','B','','',''),font.label = list(size = 32, color = "black"))  + theme(plot.margin = margin(b=0.5,l=0.1,t=0.1,r=0.1, "cm")) #incr marrgin at bottom
dev.off()


######################################
#READ COUNT PER PATIENT BY SEQUENCING CENTER -- best
######################################

PLOT_OVERALL_READ_DEPTHS_BY_RACE_AND_CENTER <- function(MAPPED_READS_OR_ALL_READS){
  
  
  if (MAPPED_READS_OR_ALL_READS=='all_reads'){DF_TO_PLOT <- ALL_READS_FOR_PLOTTING} 
  if (MAPPED_READS_OR_ALL_READS=='mapped_reads'){DF_TO_PLOT <- MAPPED_READS_FOR_PLOTTING} 
  
  #remove KIRP b/c not significant
  DF_TO_PLOT <- subset(DF_TO_PLOT, (Cancer != 'KIRP'))
  DF_TO_PLOT <- droplevels(DF_TO_PLOT)
  
  #remove any with zero counts  
  DF_TO_PLOT <- subset(DF_TO_PLOT, read_depth != 0)
  
  
  seq_site<- DF_TO_PLOT
  
  #change variable names
  seq_site$Sequencing_Center <- gsub("_"," ",seq_site$Sequencing_Center)
  
  #replace source site levels with numbers, to de-identify
  seq_site_numbers <- data.frame(unique(seq_site$Sequencing_Center))
  seq_site_numbers$Sequencing_Center_number <- rownames(seq_site_numbers)
  names(seq_site_numbers)[1] <- 'Sequencing_Center'
  seq_site <- merge(seq_site,seq_site_numbers, by='Sequencing_Center')
  
  seq_site$Sequencing_Center_number <- paste0("Center ",seq_site$Sequencing_Center_number) #add 'Center' to center #
  
  #calcualte number of black, white samples per site
  seq_site_counts <- seq_site[,c('submitter_id','race_PCA','Sequencing_Center_number','read_depth','Cancer','tissue.definition')] %>%
    #  tidyr::complete(sample_name, HLA_type,peptide_class, fill = list(number_of_unique_alleles = 0)) %>% #retain zero-counts
    dplyr::group_by(Cancer,tissue.definition,Sequencing_Center_number,race_PCA,.drop=FALSE) %>%
    dplyr::summarise(number_of_samples=dplyr::n()) 
  
  seq_site_counts <- seq_site_counts %>% spread(race_PCA, number_of_samples)
  seq_site_counts <- subset(seq_site_counts, (Black >= 5 & White >= 5))
  
  seq_site <- seq_site[paste(seq_site$Cancer, seq_site$Sequencing_Center_number, seq_site$tissue.definition) %in% paste(seq_site_counts$Cancer, seq_site_counts$Sequencing_Center_number, seq_site_counts$tissue.definition),]
  
  seq_site$Race_tissue <- ifelse(seq_site$race_PCA=='Black' & seq_site$tissue.definition=='Primary Tumor','Black tumor',ifelse(seq_site$race_PCA=='White' & seq_site$tissue.definition=='Primary Tumor','White tumor',ifelse(seq_site$race_PCA=='Black' & seq_site$tissue.definition=='Patient Germline','Black germline',ifelse(seq_site$race_PCA=='White' & seq_site$tissue.definition=='Patient Germline','White germline',NA))))
  
  seq_site$Race_tissue <- factor(seq_site$Race_tissue, levels=c('White tumor','Black tumor','White germline','Black germline'),ordered=TRUE)
  
  
  
  
  
  DF_TO_PLOT <- seq_site
  
  for (NUMBER in c(names(table(DF_TO_PLOT$Sequencing_Center_number)[1]),names(table(DF_TO_PLOT$Sequencing_Center_number)[2]),
                   names(table(DF_TO_PLOT$Sequencing_Center_number)[3]))){  
    
    SEQSITES_site <- subset(DF_TO_PLOT, (Sequencing_Center_number==NUMBER))
    SEQSITES_site <- droplevels(SEQSITES_site)
    # SEQSITES_cancer$SEQSITES_number <- (factor(SEQSITES_cancer$CAPTURE_KIT_number, levels=sort(unique(as.numeric(SEQSITES_cancer$CAPTURE_KIT_number))),ordered=TRUE))
    #  SEQSITES_cancer$SEQSITES_number_Kit <- SEQSITES_cancer$Sequencing_Center_number #add 'kit' to kit #
    #  SEQSITES_cancer$SEQSITES_number_Kit <- reorder(SEQSITES_cancer$SEQSITES_number_Kit, as.numeric(SEQSITES_cancer$SEQSITES_number))#now, must reorder Center by center number!!
    
    #  if(length(unique(SEQSITES_cancer$CAPTURE_KIT_simp))==10){SCALE=1.50}
    # if(length(unique(SEQSITES_cancer$CAPTURE_KIT_simp))==9){SCALE=6.10}
    #  if(length(unique(SEQSITES_cancer$CAPTURE_KIT_simp))==8){SCALE=1}
    #  if(length(unique(SEQSITES_cancer$CAPTURE_KIT_simp))==7){SCALE=04.65}
    #  if(length(unique(SEQSITES_cancer$CAPTURE_KIT_simp))==6){SCALE=08.80}
    #  if(length(unique(SEQSITES_cancer$CAPTURE_KIT_simp))==5){SCALE=27.75}
    if(length(unique(SEQSITES_site$Cancer))==4){SCALE=2}
    if(length(unique(SEQSITES_site$Cancer))==3){SCALE=5}
    if(length(unique(SEQSITES_site$Cancer))==2){SCALE=00}
    if(length(unique(SEQSITES_site$Cancer))==1){SCALE=54.50}
    
    Counts_SC <- aggregate(read_depth ~  Race_tissue + Cancer , SEQSITES_site, length) 
    
    SEQSITES_site_plot <- ggplot(SEQSITES_site, aes(x=factor(Race_tissue), y=read_depth,fill=Race_tissue)) + xlab("") + ylab("") + facet_grid(Sequencing_Center_number~Cancer,scales='free',drop=FALSE) + #in 'fill' argument can change factor order
      
      geom_violin() + geom_sina(size=0.4) +  stat_summary(geom='crossbar',fun='median',color='white',fatten=6,width=0.35) + guides(fill = guide_legend(override.aes = list(linetype = 0))) + 
      theme(plot.margin=unit(c(-0.5,SCALE,-0.5,0), "cm"),
            plot.title = element_text(size=46,hjust=0.5),
            strip.text.x = element_text(size=40),
            strip.text.y = element_text(size=40),
            axis.title.y = element_text(size=30,margin = margin(t = 0, r = 5, b = 0, l = 0)),
            axis.title.x = element_text(margin = margin(t = 15, r = 10, b = 5, l = 10)),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size=34,margin=margin(r=1)),
            axis.ticks.x = element_blank(),
            strip.background = element_blank(),
            panel.background = element_rect(color = 'black',fill=NA),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.spacing.x=unit(0.65,"cm"),
            legend.key.size = unit(4,"line"),
            legend.text=element_text(size=27), legend.title=element_blank(),
            legend.background = element_rect(linetype='solid', colour='black')) + guides(color = guide_legend(override.aes = list(size=30))) +
      
      
      
      
      
      #legend.position = c(0.5,-0.04),legend.direction='horizontal',plot.margin = margin(t=5,r=15,b=15,l=5))+
      #    geom_point() + #stat_summary(fun.data=n_fun,geom='text',size=11) + #calculate n
      
      stat_compare_means(size=11,method='wilcox.test',mapping=aes(label=..p.signif..),
                         comparisons = list(c("Black tumor","White tumor")),bracket.size = NA,tip.length = 0,color='black',
                         symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("p<0.0001", "p<0.001", "p<0.01", "p<0.05", ""))) +
      
      stat_compare_means(size=11,method='wilcox.test',mapping=aes(label=..p.signif..),
                         comparisons = list(c("Black germline","White germline")),bracket.size = NA,tip.length = 0,color='black',
                         symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("p<0.0001", "p<0.001", "p<0.01", "p<0.05", ""))) +
      
      
      scale_y_continuous(expand = expansion(mult = c(0.07, 0.235)),labels=unit_format(unit = "",scale = 1e-6,accuracy=1),trans=pseudo_log_trans(base=2),breaks = c(25000000,50000000,100000000,200000000,400000000)) +  scale_x_discrete(drop=FALSE) + #drop false to keep blank spots
      scale_fill_manual(values=c('#EF5350','#42A5F5','#FFCDD2','#BBDEFB'),labels=c('White tumor  ','Black tumor  ','White germline  ','Black germline   ')) + #MANUALLY SET FILL COLORS 
      
      
      
      coord_cartesian(ylim=c(15000000,NA),clip='off') + 
      
      #add text annotation
      {if (NUMBER=='Center 1'){geom_text(data = Counts_SC,aes(label=read_depth, y=min(SEQSITES_site$read_depth)-5500000),size=9) }} + 
      {if (NUMBER!='Center 1'){geom_text(data = Counts_SC,aes(label=read_depth, y=min(SEQSITES_site$read_depth)-6500000),size=9) }}  
    
    assign(paste0("plot_",NUMBER), SEQSITES_site_plot)
    
  }
  
  
  #combine plots and output
  
  figure <- ggarrange(NULL,get(paste0("plot_",names(table(DF_TO_PLOT$Sequencing_Center_number)[1]))),NULL,get(paste0("plot_",names(table(DF_TO_PLOT$Sequencing_Center_number)[2]))),NULL,get(paste0("plot_",names(table(DF_TO_PLOT$Sequencing_Center_number)[3]))),NULL,nrow=7,common.legend = TRUE, legend="bottom",heights=c(0.03,1,0.07,1,0.07,1,0.15))
  figure <- annotate_figure(figure,left = text_grob("Millions of reads per exome (log2 scale)",size=40,rot = 90)) + theme(plot.margin = margin(b=0.75,l=0.5,t=0.75,r=2, "cm")) #incr marrgin at bottom
  
  
}


SEQ_CENTER_MAPPED_READS <- PLOT_OVERALL_READ_DEPTHS_BY_RACE_AND_CENTER('mapped_reads')

#OUTPUT FIGURE S3
#pdf(paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/GitHub/Racial_Disparities_in_Exome_Read_Depth/MS_plots/FIGS3.pdf'),height=20,width=15,onefile = FALSE)
#ggarrange(SEQ_CENTER_MAPPED_READS)
#dev.off()

#OUTPUT FIGURE S3
png(paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/GitHub/Racial_Disparities_in_Exome_Read_Depth/MS_plots/FIGS3.png'),height=20,width=15,res=500,units='in')
ggarrange(SEQ_CENTER_MAPPED_READS)
dev.off()





######################################
#READ COUNT PER PATIENT BY SAMPLE COLLECTION SITE
######################################
PLOT_OVERALL_READ_DEPTHS_BY_RACE_AND_SITE <- function(MAPPED_READS_OR_ALL_READS){
  
  
  if (MAPPED_READS_OR_ALL_READS=='all_reads'){DF_TO_PLOT <- ALL_READS_FOR_PLOTTING} 
  if (MAPPED_READS_OR_ALL_READS=='mapped_reads'){DF_TO_PLOT <- MAPPED_READS_FOR_PLOTTING} 
  
  #remove KIRP, as nothing was significant
  DF_TO_PLOT <- subset(DF_TO_PLOT,(Cancer != 'KIRP'))
  DF_TO_PLOT <- droplevels(DF_TO_PLOT)
  
  #remove any with zero counts  
  DF_TO_PLOT <- subset(DF_TO_PLOT, read_depth != 0)
  
  source_site <- DF_TO_PLOT
  
  
  
  #calcualte number of black, white samples per site
  source_site_counts <- source_site[,c('submitter_id','race_PCA','Source_Site','read_depth','Cancer','tissue.definition')] %>%
    #  tidyr::complete(sample_name, HLA_type,peptide_class, fill = list(number_of_unique_alleles = 0)) %>% #retain zero-counts
    dplyr::group_by(Cancer,tissue.definition,Source_Site,race_PCA,.drop=FALSE) %>%
    dplyr::summarise(number_of_samples=dplyr::n()) 
  
  
  #'spread' the data and only keep centers with at least 5 black and at least 5 white samples
  source_site_counts <- source_site_counts %>% spread(race_PCA, number_of_samples)
  source_site_counts <- data.frame(subset(source_site_counts, (Black >= 5 & White >= 5)))
  
  #replace source site levels with numbers, to de-identify
  source_site_numbers <- data.frame(unique(source_site_counts$Source_Site))
  source_site_numbers$Source_Site_number <- rownames(source_site_numbers)
  names(source_site_numbers)[1] <- 'Source_Site'
  source_site_counts <- merge(source_site_counts,source_site_numbers, by='Source_Site')
  
  source_site <- source_site[paste(source_site$Cancer, source_site$Source_Site, source_site$tissue.definition) %in% paste(source_site_counts$Cancer, source_site_counts$Source_Site, source_site_counts$tissue.definition),]
  
  source_site <- merge(source_site,source_site_numbers, by='Source_Site')
  
  
  #new column containing race_PCA and tissue
  source_site$Race_tissue <- ifelse(source_site$race_PCA=='Black' & source_site$tissue.definition=='Primary Tumor','Black tumor',ifelse(source_site$race_PCA=='White' & source_site$tissue.definition=='Primary Tumor','White tumor',ifelse(source_site$race_PCA=='Black' & source_site$tissue.definition=='Patient Germline','Black germline',ifelse(source_site$race_PCA=='White' & source_site$tissue.definition=='Patient Germline','White germline',NA))))
  
  source_site$Race_tissue <- factor(source_site$Race_tissue, levels=c('White tumor','Black tumor','White germline','Black germline'))
  
  source_site$Source_Site_number2 <- factor(source_site$Source_Site_number,levels=paste0(seq(1:19)))
  source_site$Source_Site_number3 <- paste0("Site ",source_site$Source_Site_number2)
  source_site$Source_Site_number3 <- factor(source_site$Source_Site_number3, levels=paste0("Site ",seq(1:19)),ordered=TRUE)
  
  
  #for BRCA, create two new...
  source_site$Cancer <- ifelse((source_site$Cancer=='BRCA' & (source_site$Source_Site_number3=='Site 6' | source_site$Source_Site_number3=='Site 7' | source_site$Source_Site_number3=='Site 8' | source_site$Source_Site_number3=='Site 9')),as.character('BRCA '),as.character(source_site$Cancer))
  source_site$Cancer <- factor(source_site$Cancer, levels=c('BRCA','BRCA ','UCEC','PRAD','LUAD','COAD','KIRC'))
  
  
  
  
  for (CANCER_NAME in c('BRCA','BRCA ','UCEC','PRAD','LUAD','COAD','KIRC')){  
    
    source_site_cancer <- subset(source_site, (Cancer==CANCER_NAME))
    source_site_cancer <- droplevels(source_site_cancer)
    
    
    
    #    source_site_cancer$Source_Site_number <- reorder(source_site_cancer$Source_Site_number, as.numeric(source_site_cancer$Source_Site_number))#now, must reorder Center by center number!!
    
    
    #  if(length(unique(source_site_cancer$Source_Site_number))==9){SCALE=6.10}
    #  if(length(unique(source_site_cancer$Source_Site_number))==8){SCALE=1}
    #  if(length(unique(source_site_cancer$Source_Site_number))==7){SCALE=37}
    #  if(length(unique(source_site_cancer$Source_Site_number))==6){SCALE=08.80}
    #  if(length(unique(source_site_cancer$Source_Site_number))==5){SCALE=55.5}
    # if(length(unique(source_site_cancer$Source_Site_number))==4){SCALE=76}
    #  if(length(unique(source_site_cancer$Source_Site_number))==3){SCALE=80.0}
    #  if(length(unique(source_site_cancer$Source_Site_number))==2){SCALE=92.25}
    #  if(length(unique(source_site_cancer$Source_Site_number))==1){SCALE=104.60}
    
    if(length(unique(source_site_cancer$Source_Site_number))==5){SCALE=1}
    if(length(unique(source_site_cancer$Source_Site_number))==4){SCALE=16.5}
    if(length(unique(source_site_cancer$Source_Site_number))==3){SCALE=32.1}
    if(length(unique(source_site_cancer$Source_Site_number))==2){SCALE=47.8}
    if(length(unique(source_site_cancer$Source_Site_number))==1){SCALE=63.25}
    
    
    Counts_CS <- aggregate(read_depth ~  Source_Site_number3 + Race_tissue , source_site_cancer, length)   
    
    #get range of values for y axis
    RANGE1 <- round(range(source_site_cancer$read_depth),digits = -5 )[1]
    RANGE2 <- round(range(source_site_cancer$read_depth),digits = -5 )[2]
    
    
    source_site_cancer_plot <- ggplot(source_site_cancer, aes(x=factor(Race_tissue), y=read_depth,fill=Race_tissue)) + xlab("") + ylab("") + facet_grid(Cancer~Source_Site_number3,scales='free',drop=FALSE) + #in 'fill' argument can change factor order
      
      geom_violin() + geom_sina(size=1) +  stat_summary(geom='crossbar',fun='median',color='white',fatten=6,width=0.35) + guides(fill = guide_legend(override.aes = list(linetype = 0))) + 
      labs(title = "") +
      theme(plot.margin=unit(c(-0.5,SCALE,-0.5,0), "cm"),
            plot.title = element_text(size=46,hjust=0.5),
            strip.text.x = element_text(size=52),
            strip.text.y = element_text(size=52),
            axis.title.y = element_text(size=30,margin = margin(t = 0, r = 5, b = 0, l = 0)),
            axis.title.x = element_text(margin = margin(t = 15, r = 10, b = 5, l = 10)),
            axis.text.x = element_blank(),
            axis.ticks.x= element_blank(),
            axis.text.y = element_text(size=42,margin=margin(r=1)),
            strip.background = element_blank(),
            panel.background = element_rect(color = 'black',fill=NA),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.spacing.x=unit(0.65,"cm"),
            legend.key.size = unit(4,"line"),
            legend.text=element_text(size=70), legend.title=element_blank(),
            legend.background = element_rect(linetype='solid', colour='black')) + guides(color = guide_legend(override.aes = list(size=30))) +
      
      
      #legend.position = c(0.5,-0.04),legend.direction='horizontal',plot.margin = margin(t=5,r=15,b=15,l=5))+
      #  geom_point() +
      
      stat_compare_means(size=15,method='wilcox.test',mapping=aes(label=..p.signif..),
                         comparisons = list(c("Black tumor","White tumor")),bracket.size = NA,tip.length = 0,color='black',
                         symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("p<0.0001", "p<0.001", "p<0.01", "p<0.05", ""))) +
      
      stat_compare_means(size=15,method='wilcox.test',mapping=aes(label=..p.signif..),
                         comparisons = list(c("Black germline","White germline")),bracket.size = NA,tip.length = 0,color='black',
                         symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("p<0.0001", "p<0.001", "p<0.01", "p<0.05", ""))) +
      
      
      scale_y_continuous(expand = expansion(mult = c(0.07, 0.19)),labels=unit_format(unit = "",scale = 1e-6,accuracy=1),trans=pseudo_log_trans(base=2),breaks = round_any(round(2^seq(from=log(RANGE1,base=2),to=log(RANGE2,base=2),length.out=5),digits = -5),10000000)) + 
      
      
      
      
      {if(CANCER_NAME=='PRAD' ){coord_cartesian(ylim=c(28000000,NA),clip='off') }} + 
      {if(CANCER_NAME=='LUAD' ){coord_cartesian(ylim=c(40000000,NA),clip='off') }} + 
      {if(CANCER_NAME=='BRCA' | CANCER_NAME=='BRCA ' | CANCER_NAME=='UCEC' ){coord_cartesian(ylim=c(30000000,NA),clip='off') }} + 
      
      {if(CANCER_NAME=='KIRC' ){coord_cartesian(ylim=c(40000000,NA),clip='off') }} + 
      {if(CANCER_NAME=='COAD' ){coord_cartesian(ylim=c(16000000,NA),clip='off') }} + 
      
      
      
      
      scale_x_discrete(drop=FALSE) + #drop false to keep blank spots
      
      
      
      scale_fill_manual(values=c('#EF5350','#42A5F5','#FFCDD2','#BBDEFB'),labels=c('White tumor  ','Black tumor  ','White germline  ','Black germline   ')) + #MANUALLY SET FILL COLORS
      
      #add text annotation
      {if(CANCER_NAME=='PRAD'){geom_text(data = Counts_CS,aes(label=read_depth, y=min(source_site_cancer$read_depth)-32500000,vjust=0),size=14)}} + 
      {if(CANCER_NAME=='LUAD'){geom_text(data = Counts_CS,aes(label=read_depth, y=min(source_site_cancer$read_depth)-19500000,vjust=0),size=14)}} + 
      {if(CANCER_NAME=='BRCA'){geom_text(data = Counts_CS,aes(label=read_depth, y=min(source_site_cancer$read_depth)-32500000,vjust=0),size=14)}} + 
      {if(CANCER_NAME=='BRCA '){geom_text(data = Counts_CS,aes(label=read_depth, y=min(source_site_cancer$read_depth)-32500000,vjust=0),size=14)}} + 
      {if(CANCER_NAME=='KIRC'){geom_text(data = Counts_CS,aes(label=read_depth, y=min(source_site_cancer$read_depth)-23500000,vjust=0),size=14)}} + 
      {if(CANCER_NAME=='COAD'){geom_text(data = Counts_CS,aes(label=read_depth, y=min(source_site_cancer$read_depth)-41000000,vjust=0),size=14)}} + 
      {if(CANCER_NAME=='UCEC'){geom_text(data = Counts_CS,aes(label=read_depth, y=min(source_site_cancer$read_depth)-31000000,vjust=0),size=14)}} 
    
    
    #now convert to grob, identify which cols contain facets, and set facet width
    # CAPKITS_cancer <- ggplot_gtable(ggplot_build(CAPKITS_cancer))
    #  facet.cols <- CAPKITS_cancer$layout$l[grepl("panel", CAPKITS_cancer$layout$name)]
    # CAPKITS_cancer$widths[facet.cols] <- unit(c(5), "cm")
    # CAPKITS_cancer$heights[facet.cols] <- unit(c(2),"null")
    assign(paste0("plot_",CANCER_NAME), source_site_cancer_plot)
    
  }
  
  
  #combine plots and output
  
  figure <- ggarrange(NULL,plot_BRCA,NULL,`plot_BRCA `,NULL,plot_UCEC,NULL,plot_PRAD,NULL,plot_LUAD,NULL,plot_COAD,NULL,plot_KIRC,NULL,nrow=15,common.legend = TRUE, legend="bottom",heights=c(0.10,1.5,0.15,1.5,0.15,1.5,0.15,1.5,0.15,1.5,0.15,1.5,0.15,1.5,0.30))
  figure <- annotate_figure(figure,left = text_grob("Millions of reads per exome (log2 scale)",size=60,rot = 90)) + theme(plot.margin = margin(b=0.75,l=0.5,t=0.75,r=2, "cm")) #incr marrgin at bottom
  
  
}


COLLECTION_SITE_MAPPED_READS <- PLOT_OVERALL_READ_DEPTHS_BY_RACE_AND_SITE('mapped_reads')


#OUTPUT FIGURE S4
#pdf(paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/GitHub/Racial_Disparities_in_Exome_Read_Depth/MS_plots/FIGS4.pdf'),height=50,width=35,onefile = FALSE)
#ggarrange(COLLECTION_SITE_MAPPED_READS)
#dev.off()

#OUTPUT FIGURE S4
png(paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/GitHub/Racial_Disparities_in_Exome_Read_Depth/MS_plots/FIGS4.png'),height=50,width=35,res=500,units='in')
ggarrange(COLLECTION_SITE_MAPPED_READS)
dev.off()




######################################
#READ COUNT PER PATIENT BY CAPTURE KIT
######################################
options(scipen=999)

PLOT_OVERALL_READ_DEPTHS_BY_RACE_AND_CAPKIT <- function(MAPPED_READS_OR_ALL_READS){
  
  
  if (MAPPED_READS_OR_ALL_READS=='all_reads'){DF_TO_PLOT <- ALL_READS_FOR_PLOTTING} 
  if (MAPPED_READS_OR_ALL_READS=='mapped_reads'){DF_TO_PLOT <- MAPPED_READS_FOR_PLOTTING} 
  
  #remove KIRP, as nothing was significant
  DF_TO_PLOT <- subset(DF_TO_PLOT,(Cancer != 'KIRP'))
  DF_TO_PLOT <- droplevels(DF_TO_PLOT)
  
  #remove any with zero counts  
  DF_TO_PLOT <- subset(DF_TO_PLOT, read_depth != 0)

  
  CAPKITS <- DF_TO_PLOT
  
  #change variable names
  # CAPKITS$CAPKITS <- gsub("_"," ",CAPKITS$CAPKITS)
  
  #calcualte number of black, white samples per site
  CAPKITS_counts <- CAPKITS[,c('submitter_id','race_PCA','CAPTURE_KIT_simp','read_depth','Cancer','tissue.definition')] %>%
    #  tidyr::complete(sample_name, HLA_type,peptide_class, fill = list(number_of_unique_alleles = 0)) %>% #retain zero-counts
    dplyr::group_by(Cancer,tissue.definition,CAPTURE_KIT_simp,race_PCA,.drop=FALSE) %>%
    dplyr::summarise(number_of_samples=dplyr::n()) 
  
  #only keep kit+cancer combos with more than 5 samples
  CAPKITS_counts <- subset(CAPKITS_counts, number_of_samples >= 5)
  
  
  #replace source site levels with numbers, to de-identify
  #  CAPKITS_numbers <- data.frame(unique(CAPKITS_counts$CAPTURE_KIT_simp))
  #  CAPKITS_numbers$CAPTURE_KIT_number <- rownames(CAPKITS_numbers)
  #  names(CAPKITS_numbers)[1] <- 'CAPTURE_KIT_simp'
  #  CAPKITS_counts <- merge(CAPKITS_counts,CAPKITS_numbers, by='CAPTURE_KIT_simp')
  
  #remove any below 5
  CAPKITS <- CAPKITS[paste(CAPKITS$Cancer, CAPKITS$CAPTURE_KIT_simp, CAPKITS$tissue.definition, CAPKITS$race_PCA) %in% paste(CAPKITS_counts$Cancer, CAPKITS_counts$CAPTURE_KIT_simp, CAPKITS_counts$tissue.definition,CAPKITS_counts$race_PCA),]
  
  # CAPKITS <- merge(CAPKITS,CAPKITS_numbers, by='CAPTURE_KIT_simp')
  
  
  #new column containing race_PCA and tissue
  CAPKITS$Race_tissue <- ifelse(CAPKITS$race_PCA=='Black' & CAPKITS$tissue.definition=='Primary Tumor','Black tumor',ifelse(CAPKITS$race_PCA=='White' & CAPKITS$tissue.definition=='Primary Tumor','White tumor',ifelse(CAPKITS$race_PCA=='Black' & CAPKITS$tissue.definition=='Patient Germline','Black germline',ifelse(CAPKITS$race_PCA=='White' & CAPKITS$tissue.definition=='Patient Germline','White germline',NA))))
  
  CAPKITS$Race_tissue <- factor(CAPKITS$Race_tissue, levels=c('White tumor','Black tumor','White germline','Black germline'))
  
  for (CANCER_NAME in c(names(table(DF_TO_PLOT$Cancer)[1]),names(table(DF_TO_PLOT$Cancer)[2]),
                        names(table(DF_TO_PLOT$Cancer)[3]),names(table(DF_TO_PLOT$Cancer)[4]),
                        names(table(DF_TO_PLOT$Cancer)[5]),names(table(DF_TO_PLOT$Cancer)[6]))){  
    
    CAPKITS_cancer <- subset(CAPKITS, (Cancer==CANCER_NAME))
    CAPKITS_cancer <- droplevels(CAPKITS_cancer)
    #  CAPKITS_cancer$CAPKITS_number <- (factor(CAPKITS_cancer$CAPTURE_KIT_number, levels=sort(unique(as.numeric(CAPKITS_cancer$CAPTURE_KIT_number))),ordered=TRUE))
    #  CAPKITS_cancer$CAPKITS_number_Kit <- paste0("Kit ",CAPKITS_cancer$CAPKITS_number) #add 'kit' to kit #
    #  CAPKITS_cancer$CAPKITS_number_Kit <- reorder(CAPKITS_cancer$CAPKITS_number_Kit, as.numeric(CAPKITS_cancer$CAPKITS_number))#now, must reorder Center by center number!!
    
    #  if(length(unique(CAPKITS_cancer$CAPTURE_KIT_simp))==10){SCALE=1.50}
    # if(length(unique(CAPKITS_cancer$CAPTURE_KIT_simp))==9){SCALE=6.10}
    #  if(length(unique(CAPKITS_cancer$CAPTURE_KIT_simp))==8){SCALE=1}
    #  if(length(unique(CAPKITS_cancer$CAPTURE_KIT_simp))==7){SCALE=04.65}
    #  if(length(unique(CAPKITS_cancer$CAPTURE_KIT_simp))==6){SCALE=08.80}
    #  if(length(unique(CAPKITS_cancer$CAPTURE_KIT_simp))==5){SCALE=27.75}
    if(length(unique(CAPKITS_cancer$CAPTURE_KIT_simp))==4){SCALE=5}
    if(length(unique(CAPKITS_cancer$CAPTURE_KIT_simp))==3){SCALE=21.50}
    #  if(length(unique(CAPKITS_cancer$CAPTURE_KIT_simp))==2){SCALE=25.80}
    if(length(unique(CAPKITS_cancer$CAPTURE_KIT_simp))==1){SCALE=54.50}
    
    Counts_CK <- aggregate(read_depth ~  CAPTURE_KIT_simp + Race_tissue , CAPKITS_cancer, length)
    
    #get range of values for y axis
    RANGE1 <- round(range(CAPKITS_cancer$read_depth),digits = -5 )[1]
    RANGE2 <- round(range(CAPKITS_cancer$read_depth),digits = -5 )[2]
    
    
    
    CAPKITS_cancer_plot <- ggplot(CAPKITS_cancer, aes(x=factor(Race_tissue), y=read_depth,fill=Race_tissue)) + xlab("") + ylab("") + facet_grid(Cancer~CAPTURE_KIT_simp,scales='free',drop=FALSE,labeller = label_wrap_gen(width=20)) + #in 'fill' argument can change factor order
      
      geom_violin() + geom_sina(size=1) +  stat_summary(geom='crossbar',fun='median',color='white',fatten=7,width=0.35) + guides(fill = guide_legend(override.aes = list(linetype = 0))) + 
      labs(title = "") +
      theme(plot.margin=unit(c(-0.5,SCALE,-0.5,0), "cm"),
            plot.title = element_text(size=46,hjust=0.5),
            strip.text.x = element_text(size=40),
            strip.text.y = element_text(size=46),
            axis.title.y = element_text(size=30,margin = margin(t = 0, r = 5, b = 0, l = 0)),
            axis.title.x = element_text(margin = margin(t = 15, r = 10, b = 5, l = 10)),
            axis.text.x = element_blank(),
            axis.ticks.x= element_blank(),
            axis.text.y = element_text(size=42,margin=margin(r=1)),
            strip.background = element_blank(),
            panel.background = element_rect(color = 'black',fill=NA),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.spacing.x=unit(0.65,"cm"),
            legend.key.size = unit(4,"line"),
            legend.text=element_text(size=62), legend.title=element_blank(),
            legend.background = element_rect(linetype='solid', colour='black')) + guides(color = guide_legend(override.aes = list(size=30))) +
      
      
      
      #legend.position = c(0.5,-0.04),legend.direction='horizontal',plot.margin = margin(t=5,r=15,b=15,l=5))+
      #  geom_point() + 
      
      stat_compare_means(size=13,method='wilcox.test',mapping=aes(label=..p.signif..),
                         comparisons = list(c("Black tumor","White tumor")),bracket.size = NA,tip.length = 0,color='black',
                         symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("p<0.0001", "p<0.001", "p<0.01", "p<0.05", ""))) +
      
      stat_compare_means(size=13,method='wilcox.test',mapping=aes(label=..p.signif..),
                         comparisons = list(c("Black germline","White germline")),bracket.size = NA,tip.length = 0,color='black',
                         symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("p<0.0001", "p<0.001", "p<0.01", "p<0.05", ""))) +
      
      
      scale_fill_manual(values=c('#EF5350','#42A5F5','#FFCDD2','#BBDEFB'),labels=c('White tumor  ','Black tumor  ','White germline  ','Black germline   ')) + #MANUALLY SET FILL COLORS
      
      #drop false to keep blank spots
      {if (CANCER_NAME != 'BRCA'){coord_cartesian(ylim=c(30000000,NA),clip='off')}} + 
      
      
      {if (CANCER_NAME == 'COAD'){coord_cartesian(ylim=c(12000000,NA),clip='off')}} + 
      
      
      {if (CANCER_NAME == 'BRCA'){coord_cartesian(ylim=c(30000000,300000000),clip='off')}} + 
      
      {if (CANCER_NAME == 'KIRC'){coord_cartesian(ylim=c(45000000,300000000),clip='off')}} + 
      
      scale_y_continuous(expand = expansion(mult = c(0.07, 0.20)),labels=unit_format(unit = "",scale = 1e-6,accuracy=1),trans=pseudo_log_trans(base=2),breaks = round_any(round(2^seq(from=log(RANGE1,base=2),to=log(RANGE2,base=2),length.out=5),digits = -5),1000000)) + 
      
      
      scale_x_discrete(drop=FALSE) + 
      
      #      {if (CANCER_NAME == 'BRCA'){coord_cartesian(ylim=c(17500000,300000000),clip='off') + scale_y_continuous(expand = expansion(mult = c(0.07, 0.25)),labels=unit_format(unit = "",scale = 1e-6,accuracy=1),trans=pseudo_log_trans(base=2),breaks = c(25000000,50000000,100000000,200000000)) +  scale_x_discrete(drop=FALSE)}} + 
      
      #add text annotation
      {if(CANCER_NAME=='BRCA'){geom_text(data = Counts_CK,aes(label=read_depth, y=min(CAPKITS_cancer$read_depth)-30000000,vjust=0),size=12)}} + 
      {if(CANCER_NAME=='COAD'){geom_text(data = Counts_CK,aes(label=read_depth, y=min(CAPKITS_cancer$read_depth)-43000000,vjust=0),size=12)}} + 
      
      {if(CANCER_NAME=='PRAD'){geom_text(data = Counts_CK,aes(label=read_depth, y=min(CAPKITS_cancer$read_depth)-30000000,vjust=0),size=12)}} + 
      
      
      {if (CANCER_NAME == 'UCEC' | CANCER_NAME == 'LUAD'){geom_text(data = Counts_CK,aes(label=read_depth, y=min(CAPKITS_cancer$read_depth)-30000000,vjust=0),size=12)}} + 
      {if (CANCER_NAME == 'KIRC' ){geom_text(data = Counts_CK,aes(label=read_depth, y=min(CAPKITS_cancer$read_depth)-19000000,vjust=0),size=12)}}
    
    
    if(length(unique(CAPKITS_cancer$CAPTURE_KIT_simp))==9){SCALE=6.10}
    
    #now convert to grob, identify which cols contain facets, and set facet width
    # CAPKITS_cancer <- ggplot_gtable(ggplot_build(CAPKITS_cancer))
    #  facet.cols <- CAPKITS_cancer$layout$l[grepl("panel", CAPKITS_cancer$layout$name)]
    # CAPKITS_cancer$widths[facet.cols] <- unit(c(5), "cm")
    # CAPKITS_cancer$heights[facet.cols] <- unit(c(2),"null")
    assign(paste0("plot_",CANCER_NAME), CAPKITS_cancer_plot)
    
  }
  
  
  #combine plots and output
  
  figure <- ggarrange(NULL,plot_UCEC,NULL,plot_PRAD,NULL,plot_LUAD,NULL,plot_COAD,NULL,plot_KIRC,NULL,nrow=11,common.legend = TRUE, legend="bottom",heights=c(0.12,1,0.12,1,0.12,1,0.12,1,0.12,1,0.20))
  figure <- annotate_figure(figure,left = text_grob("Millions of reads per exome (log2 scale)",size=60,rot = 90)) + theme(plot.margin = margin(b=0.75,l=0.5,t=0.75,r=2, "cm")) #incr marrgin at bottom
  
  
}

CAPTURE_KIT_MAPPED_READS <- PLOT_OVERALL_READ_DEPTHS_BY_RACE_AND_CAPKIT('mapped_reads')

#OUTPUT FIGURE S5
#pdf(paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/GitHub/Racial_Disparities_in_Exome_Read_Depth/MS_plots/FIGS5.pdf'),height=48,width=32,onefile = FALSE)
#ggarrange(CAPTURE_KIT_MAPPED_READS)
#dev.off()


png(paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/GitHub/Racial_Disparities_in_Exome_Read_Depth/MS_plots/FIGS5.png'),height=40,width=32,res=500,units='in')
ggarrange(CAPTURE_KIT_MAPPED_READS)
dev.off()

#####################################
#TABLE OF DATES
#####################################
DATES <- MAPPED_READS_FOR_PLOTTING[,c('Cancer','CAPTURE_KIT_simp','Sequencing_Date','race_PCA')]
DATES <- DATES %>% group_by(Cancer,CAPTURE_KIT_simp,Sequencing_Date,race_PCA) %>%count() 
DATES <- subset(DATES, Sequencing_Date != '')
DATES$Race_Count <- paste(DATES$n,DATES$race_PCA)
DATES <- DATES[,c('Cancer','CAPTURE_KIT_simp','Sequencing_Date','Race_Count')]

DATES <- aggregate(.~ Cancer+CAPTURE_KIT_simp+Sequencing_Date, DATES, toString) 

DATES$Race_Count <- gsub('White, ','White | ',DATES$Race_Count)
DATES$Race_Count <- gsub('White','W',DATES$Race_Count)
DATES$Race_Count <- gsub('Black', 'B',DATES$Race_Count)

DATES <- DATES %>% spread(Sequencing_Date,Race_Count,fill = '')

write.table(DATES,file='/research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/GitHub/Racial_Disparities_in_Exome_Read_Depth/MS_plots/seq_dates.txt',sep='\t',row.names=FALSE,quote=FALSE)



######################################
#READ COUNT PER BRCA BAM FOR NIMBLEGEN-CAPTURED SAMPLES ONLY (just dividing overall depth by capture size for each kit)
######################################

NIMBLEGEN_BY_CAPTURE_BRCA_OVERALL_READ_DEPTHS <- function(){
  
  DF_TO_PLOT <- BRCA_three_nimblegen_kits_overall_reads
  
  
  
  #Divide read depth by capture size
  DF_TO_PLOT$read_depth_per_mb_capkit <- ifelse(grepl('hg18',DF_TO_PLOT$CAPTURE_KIT_simp),DF_TO_PLOT$read_depth/35.95,ifelse(grepl('v3',DF_TO_PLOT$CAPTURE_KIT_simp),DF_TO_PLOT$read_depth/64.19,ifelse(grepl('v2',DF_TO_PLOT$CAPTURE_KIT_simp),DF_TO_PLOT$read_depth/80.05,NA))) 
  
  
  DF_TO_PLOT_cp <- DF_TO_PLOT #copy dataframe so that can have 'all kits' column
  DF_TO_PLOT_cp$CAPTURE_KIT_simp <- 'All NimbleGen kits'
  
  DF_TO_PLOT <- rbind(DF_TO_PLOT,DF_TO_PLOT_cp)
  
  DF_TO_PLOT$CAPTURE_KIT_simp <- factor(DF_TO_PLOT$CAPTURE_KIT_simp, levels=c('All NimbleGen kits','NimbleGen hg18 Exome v2','NimbleGen SeqCap EZ Exome v2','NimbleGen SeqCap EZ Exome v3'),ordered=TRUE)
  
  #calculate summary statistics
  Counts <- aggregate(read_depth_per_mb_capkit ~ tissue.definition +CAPTURE_KIT_simp + race_PCA, DF_TO_PLOT, length)
  Means <- aggregate(read_depth_per_mb_capkit ~ tissue.definition + CAPTURE_KIT_simp + race_PCA, DF_TO_PLOT, mean)
  Means$read_depth_per_mb_capkit <- sprintf("%1.2f",round(Means$read_depth_per_mb_capkit,1)/1000000)
  
  #create labels for the annotated text below each facet row
  tumor_label_counts <- data.frame(race_PCA='White',tissue.definition=c('Primary Tumor'), read_depth_per_mb_capkit=c(140000),CAPTURE_KIT_simp=c('All NimbleGen kits','NimbleGen SeqCap EZ Exome v2','NimbleGen SeqCap EZ Exome v3','NimbleGen hg18 Exome v2'),label=c('n','','',''))
  #tumor_label_counts$Cancer <- factor(tumor_label_counts$Cancer, levels=c('BRCA','UCEC','PRAD','LUAD','COAD','KIRP','KIRC'))
  germline_label_counts <- data.frame(race_PCA='White',tissue.definition=c('Patient Germline'), read_depth_per_mb_capkit=c(140000),CAPTURE_KIT_simp=c('All NimbleGen kits','NimbleGen SeqCap EZ Exome v2','NimbleGen SeqCap EZ Exome v3','NimbleGen hg18 Exome v2'),label=c('n','','',''))
  #germline_label_counts$Cancer <- factor(germline_label_counts$Cancer, levels=c('BRCA','UCEC','PRAD','LUAD','COAD','KIRP','KIRC'))
  
  tumor_label_mean <- data.frame(race_PCA='White',tissue.definition=c('Primary Tumor'), read_depth_per_mb_capkit=c(140000),CAPTURE_KIT_simp=c('All NimbleGen kits','NimbleGen SeqCap EZ Exome v2','NimbleGen SeqCap EZ Exome v3','NimbleGen hg18 Exome v2'),label=c('bar(x)','','',''))
  #tumor_label_mean$Cancer <- factor(tumor_label_mean$Cancer, levels=c('BRCA','UCEC','PRAD','LUAD','COAD','KIRP','KIRC'))
  germline_label_mean <- data.frame(race_PCA='White',tissue.definition=c('Patient Germline'), read_depth_per_mb_capkit=c(140000),CAPTURE_KIT_simp=c('All NimbleGen kits','NimbleGen SeqCap EZ Exome v2','NimbleGen SeqCap EZ Exome v3','NimbleGen hg18 Exome v2'),label=c('bar(x)','','',''))
  #germline_label_mean$Cancer <- factor(germline_label_mean$Cancer, levels=c('BRCA','UCEC','PRAD','LUAD','COAD','KIRP','KIRC'))
  
  
  
  COORD_TOP <- 155000 
  COORD_BOT <- 117500
  
  plot <- ggplot(DF_TO_PLOT, aes(x=factor(race_PCA),y=read_depth_per_mb_capkit,fill=race_PCA)) + xlab("") + ylab(paste0("Millions of reads per exome per Mb of capture region (log2 scale) ")) + facet_grid(tissue.definition~CAPTURE_KIT_simp,labeller = label_wrap_gen(width=19.5)) + 
    geom_violin() + geom_sina(size=0.4) +  stat_summary(geom='crossbar',fun='median',color='white',fatten=2.5,width=0.35) + guides(fill = guide_legend(override.aes = list(linetype = 0))) + 
    theme(strip.text.x = element_text(size = 20),#,hjust=0),
          strip.text.y = element_text(size = 22),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size=20,margin = margin(r=-1)),
          axis.title.y = element_text(size=24,margin=margin(t=0,r=20,b=0,l=0)),
          axis.ticks.x= element_blank(), 
          panel.spacing.y=unit(4.25,"lines"), #to control distance betweeen facet rows
          #     plot.tag = element_text(size = 30,face='bold'),
          #    plot.tag.position = c(-0.01, 0.98),
          strip.background = element_blank(),
          panel.background = element_rect(color = 'black',fill=NA),
          plot.title = element_text(hjust=0.5,size=28,face='bold'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.direction='horizontal',
          legend.text=element_text(size=26),
          legend.background = element_rect(linetype='solid', colour='black'),
          legend.position=c(0,-40),
          legend.title=element_blank()) + 
    # geom_point() + 
    theme(plot.margin=unit(c(0.25,0.25,1,0.25),"cm")) +
    
    scale_y_continuous(expand = expansion(mult = c(0, 0.10)),labels=unit_format(unit = "",scale = 1e-6,accuracy=.2),trans=pseudo_log_trans(base=2),breaks = c(250000,500000,1000000,2000000,4000000)) +
    
    
    
    #add text annotations
    geom_text(data = subset(Counts, tissue.definition == "Primary Tumor"),aes(label=read_depth_per_mb_capkit, y=COORD_BOT,vjust=0.1),size=6) + 
    geom_text(data = subset(Means, tissue.definition == "Primary Tumor"),aes(label=read_depth_per_mb_capkit, y=COORD_TOP,vjust=0.1),size=6) +   
    geom_text(data = subset(Counts, tissue.definition == "Patient Germline"),aes(label=read_depth_per_mb_capkit, y=COORD_BOT,vjust=0.1),size=6) + 
    geom_text(data = subset(Means, tissue.definition == "Patient Germline"),aes(label=read_depth_per_mb_capkit, y=COORD_TOP,vjust=0.1),size=6) +   
    geom_text(data= tumor_label_counts,aes(label=label,y=COORD_BOT,vjust=0.1,hjust=4.4,fontface="bold"),size=6) +
    geom_text(data= germline_label_counts,aes(label=label,y=COORD_BOT,vjust=0.1,hjust=4.4,fontface="bold"),size=6) +
    
    geom_text(data= tumor_label_mean,aes(label=label,y=COORD_TOP,vjust=0.1,hjust=5.25,fontface="bold"),size=6,parse=TRUE) +
    geom_text(data= tumor_label_mean,aes(label=label,y=COORD_TOP,vjust=0.1,hjust=5.28,fontface="bold"),size=6,parse=TRUE) +
    geom_text(data= tumor_label_mean,aes(label=label,y=COORD_TOP,vjust=0.1,hjust=5.31,fontface="bold"),size=6,parse=TRUE) +
    
    
    geom_text(data= germline_label_mean,aes(label=label,y=COORD_TOP,vjust=0.10,hjust=5.25,fontface="bold"),size=6, parse=TRUE) + 
    geom_text(data= germline_label_mean,aes(label=label,y=COORD_TOP,vjust=0.10,hjust=5.28,fontface="bold"),size=6, parse=TRUE) + 
    geom_text(data= germline_label_mean,aes(label=label,y=COORD_TOP,vjust=0.10,hjust=5.31,fontface="bold"),size=6, parse=TRUE) + 
    
    
    coord_cartesian(ylim=c(210000,6000000),clip='off') +
    
    scale_fill_manual(values=c('#EF5350','#42A5F5')) + #MANUALLY SET FILL COLOR
    
    #means comparisons
    stat_compare_means(size=7,method='wilcox.test',mapping=aes(label=..p.signif..),
                       comparisons = list(c("Black","White")),bracket.size = NA,tip.length = 0,color='black',
                       #label.y=400000000, #adjust position of p-value
                       symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("p<0.0001", "p<0.001", "p<0.01", "p<0.05", "")))
  
}

BRCA_NIMBLEGEN_CAPTURE_OVERALL_DEPTHS <- NIMBLEGEN_BY_CAPTURE_BRCA_OVERALL_READ_DEPTHS()


#OUTPUT FIGURE 2
pdf(paste0("/home/mayo/m187735/s212975.Wickland_Immunomics/GitHub/Racial_Disparities_in_Exome_Read_Depth/MS_plots/FIG2.pdf"),height=10.5,width=13,onefile = FALSE)
ggarrange(ggarrange(BRCA_NIMBLEGEN_CAPTURE_OVERALL_DEPTHS,NULL,common.legend = TRUE,legend='bottom',nrow=2,heights=c(1,0.03)), NULL,nrow=2,heights=c(1,0.02)) 
dev.off()

#OUTPUT FIGURE 2
#png(paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/GitHub/Racial_Disparities_in_Exome_Read_Depth/MS_plots/FIG2.png'),height=10.5,width=13,res=500,units='in')
#ggarrange(ggarrange(BRCA_NIMBLEGEN_CAPTURE_OVERALL_DEPTHS,NULL,common.legend = TRUE,legend='bottom',nrow=2,heights=c(1,0.03)), NULL,nrow=2,heights=c(1,0.02)) 
#dev.off()


######################################
#READ COUNT PER VARIANT POSITION (SOMATIC MAF)
######################################
#get submitter IDs and race to add to MAF
race_data <- ALL_READS_FOR_PLOTTING[c('submitter_id','race_PCA')]

#Create dataframe to hold MAF data for all cancers
MAF_all <- data.frame()

CANCER_LIST <- c('BRCA','UCEC','PRAD','LUAD','COAD','KIRC')

for (CANCER_NAME in CANCER_LIST){
  print(CANCER_NAME)
  MAF_cancer <- read.table(paste0('/research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/processing/TCGA/MAFs/somatic_hg38/',CANCER_NAME,'_somatic.maf'),header=TRUE,quote="",sep='\t')
  MAF_cancer <- MAF_cancer[,c('Tumor_Sample_Barcode','Matched_Norm_Sample_Barcode','Chromosome','Start_Position','End_Position','Hugo_Symbol','Reference_Allele','Tumor_Seq_Allele1','Tumor_Seq_Allele2','t_depth','t_ref_count','t_alt_count','n_depth','dbSNP_RS','Variant_Classification','Variant_Type','Center')]   
  MAF_cancer$submitter_id <- sub("^([^-]*-[^-]*-[^-]*).*", "\\1",MAF_cancer$Tumor_Sample_Barcode)
  #  names(MAF_cancer)[2:3] <- c('tumor_tissue_depth_MAF','normal_tissue_depth_MAF')
  MAF_cancer$Cancer <- CANCER_NAME
  write.table(file=paste0('/research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/processing/TCGA/MAFs/somatic_hg38/',CANCER_NAME,'_somatic_position_list.txt'),x=unique(MAF_cancer[,c('Chromosome','Start_Position','End_Position')]),row.names=FALSE,quote=FALSE)
  MAF_all <- rbind(MAF_all,MAF_cancer)
}

race_MAF_depths <- merge(race_data, MAF_all, by='submitter_id')
race_MAF_depths <- subset(race_MAF_depths, (race_PCA=='Black' | race_PCA=='White'))

#get rid of duplicates
race_MAF_depths <- unique(race_MAF_depths)


#melt the dataframe so that can have separate tumor and normal facets
colnames(MAF_all)[which(colnames(MAF_all)=='t_depth' | colnames(MAF_all)=='n_depth')] <- c('Primary Tumor','Patient Germline')
race_MAF_depths_melted <- reshape2::melt(race_MAF_depths, id.vars=c('submitter_id','Cancer','race_PCA','Tumor_Sample_Barcode','Matched_Norm_Sample_Barcode','Chromosome','Start_Position','End_Position','Hugo_Symbol','Reference_Allele','Tumor_Seq_Allele1','Tumor_Seq_Allele2','t_ref_count','t_alt_count','dbSNP_RS','Variant_Classification','Variant_Type','Center'))

#now interested only in tumor
race_MAF_depths_melted <- subset(race_MAF_depths_melted, variable == 't_depth')

#arrange cancer order
race_MAF_depths_melted$Cancer <- factor(race_MAF_depths_melted$Cancer,levels=c('BRCA','UCEC','PRAD','LUAD','COAD','KIRC'),ordered=TRUE)

#add an empty facet
#race_MAF_depths_melted$Facet <- ''

#plot
#function to calculate n and place below each box; include comma separator
n_fun <- function(x){return(data.frame(y=-2,label=formatC(length(x), format="f", big.mark = ",", digits=0)))}


MAF_reads_by_race <- ggplot(race_MAF_depths_melted, aes(x=factor(race_PCA),y=value,fill=race_PCA)) + xlab("") + ylab("# of reads per position") + #facet_grid(Facet~Cancer,scales='free', labeller = label_wrap_gen(width=19)) + labs(title='Coverage per somatic variant position') + 
  facet_grid(~Cancer,scales='free', labeller = label_wrap_gen(width=19)) + labs(title='Coverage per somatic variant position') + 
  geom_violin() + geom_sina(size=0.25) +  stat_summary(geom='crossbar',fun='median',color='white',fatten=2,width=0.35) + guides(fill = guide_legend(override.aes = list(linetype = 0))) + 
  theme(strip.text.x = element_text(size = 18,face='bold'),#,hjust=0),
        strip.text.y = element_text(size = 18,face='bold'),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=18),
        axis.title.y = element_text(size=20,face='bold',margin=margin(t=0,r=10,b=0,l=0)),
        axis.ticks= element_blank(), 
        plot.title = element_text(hjust=0.5,size=26,face='bold'),
        plot.tag = element_text(size = 30,face='bold'),
        plot.tag.position = c(-0.01, 0.98),
        plot.margin = unit(c(2,0,0,0), "lines"), #top, right, bottom, left
        legend.position = 'bottom',legend.direction='horizontal',
        legend.text=element_text(size=34),
        legend.background = element_rect(linetype='solid', colour='black'),
        legend.title=element_blank()) + 
  geom_point() + labs(tag = "B.") +
  stat_summary(fun.data=n_fun,geom='text',size=6) + #calculate n
  scale_fill_manual(values=c('#00BFC4','#F8766D')) + #MANUALLY SET FILL COLORS
  
  
  scale_y_continuous(trans=pseudo_log_trans(base=2),breaks = 2^seq(1,12, by = 3),
                     expand = expansion(mult = c(.050, 0.175))) +
  expand_limits(y=-2) + #expand the lower y-axis limit to make room for n=
  
  #means comparisons
  stat_compare_means(size=6,method='wilcox.test',mapping=aes(label=..p.signif..),
                     comparisons = list(c("Black","White")),bracket.size = NA,tip.length = 0,color='black',
                     #label.y=400000000, #adjust position of p-value
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("p<0.0001", "p<0.001", "p<0.01", "p<0.05", ""))) 



######################################
#MOSAIC PLOT FOR POSITIONS FROM BAM 
######################################


OUTPUT_NAME <- 'All_genes'

ALL_SAMPLES <- data.frame()

for (TISSUE in c('Primary_Tumor','Patient_Germline')){
  for (FILE in list.files(paste0('/research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/processing/TCGA/processed_data/DP_per_position_from_bams/BRCA_NEW/',TISSUE,'/',OUTPUT_NAME),pattern="*summary*",full.names = TRUE)){
    SAMPLE <- read.table(FILE,header=TRUE,sep='\t')
    ALL_SAMPLES <- rbind(ALL_SAMPLES,SAMPLE)
    
  }
}

#get total capture region size
common_capture_region <- read.table('/research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/processing/TCGA/processed_data/capture_kit/three_capture_kits_intersection_hg38_liftover.bed',header=FALSE,sep='\t')

common_capture_region_size=sum(common_capture_region[,3]-common_capture_region[,2]+1)

ALL_SAMPLES$TOTAL_POSITIONS_IN_CAPTURE_REGION <- common_capture_region_size

ALL_SAMPLES$POSITIONS_WITH_NO_COVERAGE <- ALL_SAMPLES$TOTAL_POSITIONS_IN_CAPTURE_REGION-(ALL_SAMPLES$POSITIONS_WITH_COVERAGE_1_TO_10+ALL_SAMPLES$POSITIONS_WITH_COVERAGE_11_TO_39+ALL_SAMPLES$POSITIONS_WITH_COVERAGE_40_PLUS)

ALL_SAMPLES$POSITIONS_WITH_COVERAGE_0_to_10 <- ALL_SAMPLES$POSITIONS_WITH_NO_COVERAGE + ALL_SAMPLES$POSITIONS_WITH_COVERAGE_1_TO_10


ALL_SAMPLES_MOSAIC <- ALL_SAMPLES


PLOT_MOSAIC <- function(WHICH_TISSUE){
  for (KIT in c('NimbleGen_hg18_Exome_v2','NimbleGen_SeqCap_EZ_Exome_v2','NimbleGen_SeqCap_EZ_Exome_v3')){
    
    if (KIT != 'All_kits'){ALL_SAMPLES <- subset(ALL_SAMPLES_MOSAIC,CAPTURE_KIT==KIT)}
    if (KIT == 'All_kits'){ALL_SAMPLES <- ALL_SAMPLES_MOSAIC}
    
    ALL_SAMPLES <- ALL_SAMPLES%>%
      dplyr::group_by(RACE,TISSUE) %>%
      dplyr::summarise(
        
        #DP_0_to_10 = sum(as.numeric(POSITIONS_WITH_COVERAGE_0_to_10)),
        #DP_0 = sum(as.numeric(POSITIONS_WITH_NO_COVERAGE)),
        DP_1_to_10 = sum(as.numeric(POSITIONS_WITH_COVERAGE_1_TO_10)),
        DP_11_to_39 = sum(as.numeric(POSITIONS_WITH_COVERAGE_11_TO_39)),
        DP_40_plus = sum(as.numeric(POSITIONS_WITH_COVERAGE_40_PLUS))
      ) 
    
    
    ALL_SAMPLES <- subset(ALL_SAMPLES,TISSUE==WHICH_TISSUE)
    
    ALL_SAMPLES <- ALL_SAMPLES[-c(2)]
    
    names(ALL_SAMPLES)[2:4] <- c('1-10','11-39','40+')
    
    ALL_SAMPLES_mosaic <- data.frame(t(ALL_SAMPLES[-1]))
    row.names(ALL_SAMPLES_mosaic) <- sub('DP_','',row.names(ALL_SAMPLES_mosaic))
    names(ALL_SAMPLES_mosaic) <- as.character(ALL_SAMPLES$RACE)
    
    ALL_SAMPLES_mosaic <- ALL_SAMPLES_mosaic[,c('Black','White')]
    
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
    par(mar = c(2, 0, 2, 0))  
    mosaicplot(ALL_SAMPLES_mosaic,shade=TRUE,cex.axis = 1.75,las=3,main = "")
    
    grid.echo()
    PLOT <- grid.grab()
    
    grid.draw(PLOT)
    
    grid.lines(x = unit(c(1.007, 1.007), "npc"),
               y = unit(c(0.62, 0.94), "npc"),
               gp = gpar(fill="blue",col='blue'))
    grid.lines(x = unit(c(1.009, 1.009), "npc"),
               y = unit(c(0.62, 0.94), "npc"),
               gp = gpar(fill="blue",col='blue'))
    grid.lines(x = unit(c(1.0085, 1.0085), "npc"),
               y = unit(c(0.62, 0.94), "npc"),
               gp = gpar(fill="blue",col='blue'),
               arrow = arrow(length = unit(0.30, "inches"), 
                             ends="last", type="closed"))
    grid.text("Over-representation", x = 1.065, y =.795,rot = 90,gp = gpar(col = 'blue',fontface='bold',fontsize=22))
    
    
    grid.lines(x = unit(c(1.007, 1.007), "npc"),
               y = unit(c(0.53, 0.23), "npc"),
               gp = gpar(fill="red",col='red'))
    grid.lines(x = unit(c(1.009, 1.009), "npc"),
               y = unit(c(0.53, 0.23), "npc"),
               gp = gpar(fill="red",col='red'))
    grid.lines(x = unit(c(1.0085, 1.0085), "npc"),
               y = unit(c(0.53, 0.23), "npc"),
               gp = gpar(fill="red",col='red'),
               arrow = arrow(length = unit(0.30, "inches"), 
                             ends="last", type="closed"))
    
    grid.text("Under-representation", x = 1.065, y =.365,rot = 90,gp = gpar(col = 'red',fontface='bold',fontsize=22))
    
    grid.rect(x = unit(0.94, "npc"), y = unit(0.10, "npc"),
              width = unit(0.12, "npc"), height = unit(0.2, "npc"),
              gp=gpar(col="white",fill="white")) #add a white rectangle to block out 'standardized residuals' text
    
    
    PLOT <- grid.grab()
    
    assign(paste0('plot_',KIT),PLOT)
    assign(paste0('CHISQ_',KIT,'_',WHICH_TISSUE),chisq.test(ALL_SAMPLES_mosaic),envir = .GlobalEnv)
    
  }
  MOSAIC_PLOT <- cowplot::plot_grid(NULL,plot_NimbleGen_hg18_Exome_v2,plot_NimbleGen_SeqCap_EZ_Exome_v2,plot_NimbleGen_SeqCap_EZ_Exome_v3,nrow=1,rel_widths =c(0.05,0.95,0.965,0.95),labels=c('','       NimbleGen hg18 Exome v2','NimbleGen SeqCap EZ Exome v2','   NimbleGen SeqCap EZ Exome v3'),vjust=1.4,hjust=0.055,label_size = 27)
  return(MOSAIC_PLOT)
  
}

MOSAIC_GERMLINE <- PLOT_MOSAIC('Patient_Germline')
MOSAIC_TUMOR <- PLOT_MOSAIC('Primary_Tumor')





#######################
#CHI SQUARE ANALYSIS
#######################


CHISQ_TABLE <- function(TISSUE){
  
  for (KIT in c('NimbleGen_hg18_Exome_v2','NimbleGen_SeqCap_EZ_Exome_v2','NimbleGen_SeqCap_EZ_Exome_v3')){
    CHISQ <- get(paste0('CHISQ_',KIT,'_',TISSUE))
    
    
    CHISQ_T <- matrix(mapply(function(x, y){paste(x,y, sep = " | ")}, formatC(CHISQ$observed,format='f',big.mark=',',digits=0),formatC(CHISQ$expected,format='f',big.mark=',',digits = 0)),nrow=3,ncol=2,dimnames = list(c('1-10','11-39','40+'),c('Black','White')))
    
    KIT_orig <- KIT
    KIT <- gsub("_",' ',KIT)
    
    assign(paste0('CHISQ_TABLE_',KIT_orig,'_',TISSUE),ggtexttable(CHISQ_T,
                                                                  theme= ttheme(
                                                                    base_style = "default",
                                                                    base_size = 8,
                                                                    base_colour = "black",
                                                                    padding = unit(c(0.40, 0.25), "cm"),
                                                                    colnames.style = colnames_style(size = 13,fill=c('grey70','grey70','grey70','grey70','grey70','grey70','white')),
                                                                    rownames.style = rownames_style(size = 12,face = 'plain',hjust=1), #add x= if want to change position
                                                                    tbody.style = tbody_style(size = 11) )) %>% tab_add_title(KIT,hjust = -0.02,size=14) %>%
             tab_add_footnote(
               text=paste0("X\U00B2 = ",round(CHISQ$statistic,digits=0), ', p = ',CHISQ$p.value),
               face = 'italic',
               size = NULL,
               color = NULL,
               family = NULL,
               padding = unit(1.5, "line"),
               just = "right",
               hjust = NULL,
               vjust = NULL
               
               
             ),envir = .GlobalEnv)
  }
}

CHISQ_TABLE('Patient_Germline')
CHISQ_TABLE('Primary_Tumor')


#OUTPUT FIGURE S6
png(paste0("/home/mayo/m187735/s212975.Wickland_Immunomics/GitHub/Racial_Disparities_in_Exome_Read_Depth/MS_plots/FIGS6.png"),height=11,width=12.5,units = 'in',res=500)
LEFT <- ggarrange(CHISQ_TABLE_NimbleGen_hg18_Exome_v2_Patient_Germline,CHISQ_TABLE_NimbleGen_SeqCap_EZ_Exome_v2_Patient_Germline,CHISQ_TABLE_NimbleGen_SeqCap_EZ_Exome_v3_Patient_Germline,nrow=3,common.legend = TRUE,legend = 'none',align='hv',heights=c(0.8,0.055),labels=c('',''),font.label = list(size = 22, color = "black"))
LEFT <- annotate_figure(LEFT,top = text_grob("Patient Germline",size = 16,face='bold'))

RIGHT <- ggarrange(CHISQ_TABLE_NimbleGen_hg18_Exome_v2_Primary_Tumor,CHISQ_TABLE_NimbleGen_SeqCap_EZ_Exome_v2_Primary_Tumor,CHISQ_TABLE_NimbleGen_SeqCap_EZ_Exome_v3_Primary_Tumor,nrow=3,common.legend = TRUE,legend = 'none',align='hv',heights=c(0.8,0.055),labels=c('',''),font.label = list(size = 22, color = "black"))
RIGHT <- annotate_figure(RIGHT,top = text_grob("Primary Tumor",size = 16,face='bold'))


ggarrange(LEFT,NULL,RIGHT,nrow=1,ncol=3,widths = c(1,-0.20,1)) + theme(plot.margin = margin(b=0,l=0,t=0,r=0, "cm")) #incr marrgin at bottom
dev.off()



######################################
#MINOR ALLELE FREQUENCY IN GERMLINE VS. DATABASE
######################################

###Load data containing calculated MAF and DP at each locus by different coverage levels

COMBINE_DATA <- function(PREFIX,CANCER_NAME,RACE){  
  
  MAFs_DPs_TABLE <- data.frame()  
  for (TISSUE in c("Blood_Normal")){
    
    
    print(paste0("Processing ",TISSUE,": ",RACE)) 
    
    for (i in seq(1,22)){   
      
      print(paste0("  chr ",i))  
      TABLE <- read.table(paste0('/research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/processing/TCGA/processed_data/annotate_PASS_TCGA/',PREFIX,'_all_chr/',CANCER_NAME,'/',TISSUE,'/chr',i,'/',TISSUE,'_',RACE,'_final_VQSR_PASS_population_MAFs_cohort_DP_and_MAF_by_coverage_level.txt'),header=TRUE)
      
      TABLE$Race <- RACE
      
      TABLE$Tissue <- TISSUE
      
      MAFs_DPs_TABLE <- bind_rows(MAFs_DPs_TABLE,TABLE)
      
    }
    
  }
  return(MAFs_DPs_TABLE)
}

Black_MAFs_DPs_at_ALT_positions <- COMBINE_DATA('nimblegen_capture_NOV19','BRCA','BlackWhite_Black') 
Black_MAFs_DPs_at_ALT_positions <- Black_MAFs_DPs_at_ALT_positions[-c(which(names(Black_MAFs_DPs_at_ALT_positions) == 'EUR_1kGenomes'))]
names(Black_MAFs_DPs_at_ALT_positions)[which(names(Black_MAFs_DPs_at_ALT_positions) == 'AFR_1kGenomes')] <- 'BAF_1kGenomes'

White_MAFs_DPs_at_ALT_positions <- COMBINE_DATA('nimblegen_capture_NOV19','BRCA','BlackWhite_White')   #IF DESIRE DOWNSAMPLING:BlackWhite_White_DOWNSAMPLED
White_MAFs_DPs_at_ALT_positions <- White_MAFs_DPs_at_ALT_positions[-c(which(names(White_MAFs_DPs_at_ALT_positions) == 'AFR_1kGenomes'))]
names(White_MAFs_DPs_at_ALT_positions)[which(names(White_MAFs_DPs_at_ALT_positions) == 'EUR_1kGenomes')] <- 'BAF_1kGenomes'


#isolate the non-patient columns; Also, bind the two DFs together
MAFs_DPs_at_ALT_positions <- rbind(Black_MAFs_DPs_at_ALT_positions[!(grepl('TCGA',names(Black_MAFs_DPs_at_ALT_positions)) | grepl('GT',names(Black_MAFs_DPs_at_ALT_positions)))],White_MAFs_DPs_at_ALT_positions[!(grepl('TCGA',names(White_MAFs_DPs_at_ALT_positions)) | grepl('GT',names(White_MAFs_DPs_at_ALT_positions)))])

#keep only the rows that have at least 25 patients in any of the categories
MAFs_DPs_at_ALT_positions <- subset(MAFs_DPs_at_ALT_positions, (count_patients_with_DP_below_5_5 >= 25 | count_patients_with_DP_below_10_10 >= 25 | count_patients_with_DP_within_11_39 >= 25 | count_patients_with_DP_above_40_40 >= 25))


#Now, melt!
MAFs_DPs_at_ALT_positions_MELTED <- melt(MAFs_DPs_at_ALT_positions, id.vars = names(MAFs_DPs_at_ALT_positions)[!(grepl('BAF',names(MAFs_DPs_at_ALT_positions)) & grepl('patients',names(MAFs_DPs_at_ALT_positions)))])

names(MAFs_DPs_at_ALT_positions_MELTED)[names(MAFs_DPs_at_ALT_positions_MELTED)=='variable'] <- 'DP_group'
names(MAFs_DPs_at_ALT_positions_MELTED)[names(MAFs_DPs_at_ALT_positions_MELTED)=='value'] <- 'BAF'
MAFs_DPs_at_ALT_positions_MELTED$DP_group <- gsub('BAF_','',MAFs_DPs_at_ALT_positions_MELTED$DP_group)


MAFs_DPs_at_ALT_positions_MELTED <- melt(MAFs_DPs_at_ALT_positions_MELTED, id.vars = names(MAFs_DPs_at_ALT_positions_MELTED)[!(grepl('count',names(MAFs_DPs_at_ALT_positions_MELTED)) & grepl('patients',names(MAFs_DPs_at_ALT_positions_MELTED)))])

names(MAFs_DPs_at_ALT_positions_MELTED)[names(MAFs_DPs_at_ALT_positions_MELTED)=='variable']  <- 'DP_group2'
names(MAFs_DPs_at_ALT_positions_MELTED)[names(MAFs_DPs_at_ALT_positions_MELTED)=='value']  <- 'count'
MAFs_DPs_at_ALT_positions_MELTED$DP_group2 <- gsub('count_','',MAFs_DPs_at_ALT_positions_MELTED$DP_group2)


#only want the BAFs that include the ./. genotypes as 0...
MAFs_DPs_at_ALT_positions_MELTED <- MAFs_DPs_at_ALT_positions_MELTED[!grepl("patients_with",MAFs_DPs_at_ALT_positions_MELTED$DP_group),]
MAFs_DPs_at_ALT_positions_MELTED$DP_group<- gsub("_incl_..","",MAFs_DPs_at_ALT_positions_MELTED$DP_group)


#BAFs and counts should match up by group
MAFs_DPs_at_ALT_positions_MELTED <- subset(MAFs_DPs_at_ALT_positions_MELTED, DP_group==DP_group2)

#remove any below 50
MAFs_DPs_at_ALT_positions_MELTED <- subset(MAFs_DPs_at_ALT_positions_MELTED, count >= 25)  
MAFs_DPs_at_ALT_positions_MELTED <- MAFs_DPs_at_ALT_positions_MELTED[,names(MAFs_DPs_at_ALT_positions_MELTED)!='DP_group2']


#consider only biallelic, b/c multiallelic not counted properly
#MAFs_DPs_at_ALT_positions_MELTED <- subset(MAFs_DPs_at_ALT_positions_MELTED, (value <= 1))

#reorder
DP_values <- as.character(unique(MAFs_DPs_at_ALT_positions_MELTED$DP_group))
DP_values <- c(DP_values[1],rev(DP_values)[1:(length(DP_values)-1)])


MAFs_DPs_at_ALT_positions_MELTED$DP_group <- factor(MAFs_DPs_at_ALT_positions_MELTED$DP_group, levels=DP_values,ordered=TRUE)

#cosmetic alterations
MAFs_DPs_at_ALT_positions_MELTED$Race <- gsub("BlackWhite_","",MAFs_DPs_at_ALT_positions_MELTED$Race)


#also get cancer genes
#list of cancer genes from onkoKB
cancer_genes <- data.frame(read.table('/home/mayo/m187735/s212975.Wickland_Immunomics/TCGA_metadata/onkoKB_cancer_gene_list.tsv',sep='\t',header=FALSE,skip=1)[,c(1)])
names(cancer_genes) <- c('name')

#subset
MAFs_DPs_at_ALT_positions_MELTED_cancer_genes <- MAFs_DPs_at_ALT_positions_MELTED[MAFs_DPs_at_ALT_positions_MELTED$Gene %in% cancer_genes$name,] 




#plotting function
plot_MAF_by_coverage_level <- function(DF,INCLUDE_EQNS){
  
  DF$DP_group <- gsub("DP","coverage",DF$DP_group)
  DF$DP_group <- gsub("above_40_40","≥_40",DF$DP_group)
  DF$DP_group <- gsub("below_5_5","≤_5",DF$DP_group)
  DF$DP_group <- gsub("below_10_10","≤_10",DF$DP_group)
  DF$DP_group <- gsub("within_11_39","11-39",DF$DP_group)
  
  DF$Race <- gsub('BlackWhite_','',DF$Race)
  DF$Race <- gsub('_DOWNSAMPLED','',DF$Race)
  DF$Race <- factor(DF$Race, levels=c('White','Black'),ordered=TRUE)
  
  
  #PLOT
  #colors <- c('orange','limegreen','deeppink','dodgerblue','brown','darkorchid2')
  #colors <- c('#E41A1C','#377EB8','#FF7F00','#984EA3','grey20','#4DAF4A')
  colors <- c('#377EB8','grey20','#FF7F00','#984EA3')#,'#4DAF4A')
  
  #MAFs_DPs_at_ALT_positions_MELTED <- subset(MAFs_DPs_at_ALT_positions_MELTED, Gene=='TP53')) #subset!
  
  
  DF$variable2 <- ifelse(as.character(DF$DP_group) == 'all_patients',1,ifelse(as.character(DF$DP_group) == 'patients_with_coverage_≥_40',2,ifelse(as.character(DF$DP_group) == 'patients_with_coverage_11-39',3,ifelse(as.character(DF$DP_group) == 'patients_with_coverage_≤_10',4,ifelse(as.character(DF$DP_group) == 'patients_with_coverage_≤_5',5,NA)))))
  #unique(DF2$variable2)
  
  DF <- DF[order(DF$variable2) , ]
  
  DF <- subset(DF, DP_group != 'all_patients')
  
  DF$Tissue <- sub("Primary_Tumor","Primary Tumor",DF$Tissue)
  DF$Tissue <- sub("Blood_Normal","Patient Germline",DF$Tissue)
  DF$DP_group <- gsub("BAF_","",DF$DP_group)
  DF$DP_group <- gsub("_"," ",DF$DP_group)
  DF$DP_group <- gsub("patients with "," ",DF$DP_group)
  
  DP_values <- gsub("BAF_","",DP_values)
  DP_values <- gsub("_"," ",DP_values)
  DP_values <- gsub("above",">",DP_values)
  DP_values <- gsub("below","<",DP_values)
  DP_values <- gsub("patients with DP"," coverage",DP_values)
  
  #germline only
  DF <- subset(DF, Tissue == 'Patient Germline')
  
  #DF <- subset(MAFs_DPs_at_ALT_positions_MELTED, (as.numeric(as.character(value)) > 0.3 & as.numeric(as.character(BAF_1kGenomes) > 0.3)))
  #DF <- subset(MAFs_DPs_at_ALT_positions_MELTED, Gene=='BRCA2')
  
  DF$DP_group <- as.factor(as.character(fct_reorder(DF$DP_group, DF$variable2, min)))
  
  
  plot <- ggplot(DF, aes(x=as.numeric(as.character(BAF_1kGenomes)), y=as.numeric(as.character(BAF)))) + #facet_grid(rows=vars(Tissue),cols = vars(Race)) + 
    #  geom_smooth(aes(group=race,color=race),method='loess',se=FALSE,size=0.65) +
    #  facet_wrap(~ Race) +
    geom_point(shape=1,size=0.5,alpha=0.05,aes(color=(fct_reorder(DF$DP_group, DF$variable2, min)))) + ylab('TCGA BRCA cohort MAF') + xlab("1000 Genomes population MAF") + labs(title="Minor allele frequency in common capture region:\nTCGA BRCA germline vs. 1000 Genomes") + guides(color=guide_legend(title="Coverage",size=12)) + 
    
    #  geom_smooth(data=subset(DF,Race=='White'),linetype='solid',position=position_jitter(0.000),method='lm', formula= y~x,aes(color=(fct_reorder(subset(DF, Race=='White')$DP_group, subset(DF, Race=='White')$variable2, min))),se=FALSE,size=1.5)   +
    # geom_smooth(data=subset(DF,Race=='Black'),linetype='dashed',position=position_jitter(0.000),method='lm', formula= y~x,aes(color=(fct_reorder(subset(DF, Race=='Black')$DP_group, subset(DF, Race=='Black')$variable2, min))),se=FALSE,size=1.5)   +
    geom_smooth(data=DF,position=position_jitter(0.000),method='lm', formula= y~x,aes(linetype=Race,color=(fct_reorder(DF$DP_group, DF$variable2, min))),se=FALSE,size=1.5)   +
    
    #correlation equations
    #     stat_cor(size=4,aes(color=as.character(DP_group),label = paste(sub("R = ","R=",..r.label..))),  digits=4,output.type="text",  label.y = 1.15, data = subset(DF, (DP_group==DP_values[2]))) +
    #      stat_cor(size=4,aes(color=DP_group,label = paste(sub("R = ","R=",..r.label..))),  digits=4,output.type="text",  label.y = 1.12,data = subset(DF, (DP_group==DP_values[3]))) +
    #    stat_cor(size=4,aes(color=DP_group,label = paste(sub("R = ","R=",..r.label..))),  digits=4,output.type="text",  label.y = 1.09,data = subset(DF, (DP_group==DP_values[4]))) +
    #    stat_cor(size=4,aes(color=DP_group,label = paste(sub("R = ","R=",..r.label..))),  digits=4,output.type="text",  label.y = 1.06,data = subset(DF, (DP_group==DP_values[5]))) +
    #    stat_cor(size=4,aes(color=DP_group,label = paste(sub("R = ","R=",..r.label..))),  digits=4,output.type="text",  label.y = 1.03,data = subset(DF, (DP_group==DP_values[6]))) 
    
    theme(strip.text.x = element_text(size = 30),#,hjust=0),
          strip.text.y = element_text(size = 20),
          axis.text.x = element_text(size=27,margin=margin(t=5)),
          axis.text.y = element_text(size=27,margin=margin(r=5)),
          axis.title.y = element_text(size=30,margin=margin(t=0,r=20,b=0,l=0)),
          axis.title.x = element_text(size=30,margin=margin(t=20,r=0,b=0,l=0)),
          #      axis.ticks.x= element_blank(), 
          axis.ticks.length = unit(.25, "cm"),
          plot.title = element_text(hjust=0.5,size=28,face='bold'),
          plot.margin = unit(c(0,1,0,0), "lines"), #top, right, bottom, left
          plot.tag = element_text(size = 22,face='bold'),
          plot.tag.position = c(0.01, 0.98),
          legend.position=c(0.11,0.75),legend.direction='vertical',
          legend.text=element_text(size=25,hjust = 0),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.key.width =  unit(2,"line"),
          legend.key.height = unit(1,"line"),
          legend.background = element_rect(linetype='solid', colour='black'),
          panel.background = element_rect(color = 'black',fill=NA),
          legend.title=element_text(size=24)) + #labs(tag = "B.") +
    
    guides(colour = guide_legend(override.aes = list(size=2))) + #make legend lines thicker 
    guides(linetype =  guide_legend(override.aes = list(color = 'black'))) + #black color for linetype legend
    
    # ylim(-50000000,450000000) + 
    scale_color_manual(values=colors) + #specify fill colors
    
    scale_colour_manual(name = "Coverage", values = colors, 
                        labels = expression("" >= 40, "11-39", "" <= 10,"" <= 5)) + 
    #scale_linetype_manual(name = "Coverage", values = colors, 
    #                    labels = c("solid","dashed","solid","dashed")) + 
    
    
    scale_y_continuous(limits = c(0,1.0),breaks = seq(0,1, by = 0.2),expand = expansion(mult = c(0,0))) +
    scale_x_continuous(breaks = seq(0,1,by=0.2),expand = expansion(mult = c(0,0))) 
  
  {if (INCLUDE_EQNS == 'Yes'){ #if condniiton for ggplot
    #regression equations    
    plot <- plot + 
      stat_regline_equation(color=colors[1],size=4,label.y = 0.30, label.x=0.60, data = subset(DF, (DP_group==' coverage ≥ 40'))) +
      stat_regline_equation(color=colors[2],size=4,label.y = 0.25, label.x=0.60, data = subset(DF, (DP_group==' coverage 11-39'))) +
      stat_regline_equation(color=colors[3],size=4,label.y = 0.20, label.x=0.60, data = subset(DF, (DP_group==' coverage ≤ 10'))) +
      stat_regline_equation(color=colors[4],size=4,label.y = 0.15, label.x=0.60, data = subset(DF, (DP_group==' coverage ≤ 5'))) #+
    #  stat_regline_equation(color=colors[5],size=4,label.y = 0.10, label.x=0.60, data = subset(DF, (DP_group==DP_values[6])))#  stat_regline_equation(color=colors[5],size=4,label.y = 0.10, label.x=0.60, data = subset(DF, (DP_group==DP_values[6])))#  stat_regline_equation(color=colors[5],size=4,label.y = 0.10, label.x=0.60, data = subset(DF, (DP_group==DP_values[6])))#  stat_regline_equation(color=colors[5],size=4,label.y = 0.10, label.x=0.60, data = subset(DF, (DP_group==DP_values[6])))#  stat_regline_equation(color=colors[5],size=4,label.y = 0.10, label.x=0.60, data = subset(DF, (DP_group==DP_values[6])))
  }}
  return(plot)
  
  
  
}

plot.new()
LINE_PLOT_WITHOUT_EQNS <- plot_MAF_by_coverage_level(MAFs_DPs_at_ALT_positions_MELTED,'No')
#plot.new()
#LINE_PLOT_WITH_EQNS <- plot_MAF_by_coverage_level(MAFs_DPs_at_ALT_positions_MELTED,'Yes')








######################################
#THOUSAND GENOMES ALLELE FREQS WITHIN CAPTURE REGION
######################################

#load data
Thousand_Genomes_Common_capture <- read.table('/research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/processing/TCGA/processed_data/capture_kit/1000Genomes_within_common_capture_region/1000Genomes_within_common_capture_region_all_chr.vcf',sep='\t',header = FALSE)[,1:5]

#extract EUR and AFR allele frequencies
Thousand_Genomes_Common_capture$European <- str_extract(Thousand_Genomes_Common_capture$V5, "(EUR_AF=)[^;]*(?=;|$)")
Thousand_Genomes_Common_capture$African <- str_extract(Thousand_Genomes_Common_capture$V5, "(AFR_AF=)[^;]*(?=;|$)")

#isolate numbers
Thousand_Genomes_Common_capture$European <- gsub("EUR_AF=","",Thousand_Genomes_Common_capture$European)
Thousand_Genomes_Common_capture$African <- gsub("AFR_AF=","",Thousand_Genomes_Common_capture$African)

#identify positions polymorphic in at least one group
Thousand_Genomes_Common_capture <- subset(Thousand_Genomes_Common_capture, !(European ==0 & African == 0))

#cleanup
Thousand_Genomes_Common_capture <- Thousand_Genomes_Common_capture[-c(3,4,5)]
Thousand_Genomes_Common_capture$CHR_POS <- paste0(Thousand_Genomes_Common_capture$V1,"_",Thousand_Genomes_Common_capture$V2)

#rearrange dataframe for plotting
Thousand_Genomes_Common_capture_MELTED <- melt(Thousand_Genomes_Common_capture[-c(1,2)], id.vars=c('CHR_POS'))

Thousand_Genomes_Common_capture_MELTED$value <- as.numeric(Thousand_Genomes_Common_capture_MELTED$value)
Thousand_Genomes_Common_capture_MELTED$variable <- factor(Thousand_Genomes_Common_capture_MELTED$variable, levels=c('European','African'))


###################
##BOXPLOT FOR MEDIAN
###################

Means_TG <- aggregate(value ~ variable, Thousand_Genomes_Common_capture_MELTED, mean)

Mean_TG_label <- data.frame(label=c('bar(x)',''),variable=c('European','African'))

THOU_GENOMES_CAPTURE <- ggplot(Thousand_Genomes_Common_capture_MELTED, aes(x=as.factor(variable),y=value*100,fill=variable)) +  
  geom_boxplot(fatten=4) + ylab('MAF (log2scale)') + xlab('') + labs(title='Minor allele frequency in common capture region:\n1000 Genomes European vs. 1000 Genomes African') + #geom_point() + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=27,margin=margin(r=8)),
        axis.title.y = element_text(size=30,margin=margin(t=0,r=10,b=0,l=0)),
        axis.ticks.x = element_blank(), 
        plot.title = element_text(hjust=0.5,size=28,face='bold'),
        #  plot.tag = element_text(size = 30,face='bold'),
        # plot.tag.position = c(-0.01, 0.98),
        plot.margin = unit(c(0,0,0,2), "lines"), #top, right, bottom, left
        panel.background = element_rect(color = 'black',fill=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'bottom',legend.direction='horizontal',
        legend.box.margin= margin(t=25), #move margin up or down relative to plot
        legend.background = element_rect(linetype='solid', colour='black'),
        legend.text=element_text(size=31),        
        legend.title=element_blank()) + #labs(tag = "C.") +
  
  
  #coord_cartesian(ylim = c(0,1.35))+
  scale_fill_manual(values=c('#EF5350','#42A5F5')) + #MANUALLY SET FILL COLORS 
  
  #means comparisons
  stat_compare_means(size=10,method='wilcox.test',mapping=aes(label=..p.signif..),
                     comparisons = list(c("European","African")),bracket.size = NA,tip.length = 0,color='black',
                     label.y=6.5, #adjust position of p-value
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("p<0.0001", "p<0.001", "p<0.01", "p<0.05", ""))) + 
  
  #add text annotations
  geom_text(data = subset(Means_TG, variable == "European"),aes(label=round(mean(value),digits=3), y=-1.15,vjust=0),size=10) + 
  geom_text(data = subset(Means_TG, variable == "African"),aes(label=round(mean(value),digits=3), y=-1.15,vjust=0),size=10) + 
  geom_text(data= Mean_TG_label,aes(label=label,y=-1.15,vjust=0,hjust=8.15,fontface="bold"),size=11,parse=TRUE) + 
  
  
  
  scale_y_continuous(trans=scales::pseudo_log_trans(base=2) ,breaks=c(1,5,15,50),labels=function(x) paste0(x/100)) + #modify labels to have decimals
  expand_limits(y=c(0,0.2)) +
  
  coord_cartesian(ylim=c(0,NA),clip='off')



#OUTPUT FIGURE 3 (if no downsampling)
pdf(paste0("/home/mayo/m187735/s212975.Wickland_Immunomics/GitHub/Racial_Disparities_in_Exome_Read_Depth/MS_plots/FIG3.pdf"),height=19,width=24)
ROW_1 <- ggarrange(ggarrange(NULL,MOSAIC_GERMLINE,nrow=1,widths=c(0.02,1)),NULL,widths=c(1,0.05),font.label=list(size=32,color='black'),labels=c('A'))
ROW_3 <- ggarrange(NULL,LINE_PLOT_WITHOUT_EQNS,NULL,THOU_GENOMES_CAPTURE,NULL,widths=c(0.035,0.510,0.01,0.44,0.03),nrow=1,font.label=list(size=32,color='black'),labels=c('B','','','C',''))
ggarrange(NULL,ROW_1,NULL,ROW_3,NULL,widths=c(1,0.1,1,0.05),nrow=5,heights=c(0.01,0.45,0.025,0.55,0.01))
dev.off()


#OUTPUT FIGURE S63 (downsampling)
png(paste0("/home/mayo/m187735/s212975.Wickland_Immunomics/GitHub/Racial_Disparities_in_Exome_Read_Depth/MS_plots/FIGS6.png"),height=10,width=12,units='in',res=500)
ggarrange(NULL,LINE_PLOT_WITHOUT_EQNS,NULL,widths=c(0.03,1,0.03),nrow=1)
dev.off()


#OUTPUT FIGURE S
#png(paste0("/home/mayo/m187735/s212975.Wickland_Immunomics/GitHub/Racial_Disparities_in_Exome_Read_Depth/MS_plots/FIG3.png"),height=19,width=24,units='in',res=400)
#ggarrange(NULL,ROW_1,NULL,ROW_3,NULL,widths=c(1,0.1,1,0.05),nrow=5,heights=c(0.01,0.45,0.025,0.55,0.01))
#dev.off()




######################################
#BAC PER VARIANT POSITION FOR EACH SAMPLE (BRCA) - from MAF, protein-altering changes
######################################


GENERATE_BAC_PLOT <- function(MIN_READS,DATASET,X_POS,Y_POS){
  
  MINIMUM_READS <- MIN_READS
  
  MAF_subset <- DATASET
  MAF_subset$MIN_READS <- paste0(MINIMUM_READS,"x coverage") 
  MAF_subset$TYPE <- ''
  
  MAF_subset <- subset(MAF_subset, t_alt_count > 0)
  
  ##normalize b/c truncated distribution
  #MAF_subset$BAC <- trunc2norm(MAF_subset$BAC)
  
  #first create shading regions
  AUC <- density(MAF_subset$BAC, from = min(MAF_subset$BAC), to = 1/MINIMUM_READS,n=2500)
  AUC_coords <- data.frame(x = AUC$x, y = AUC$y)
  
  #create coordinates for annotation
  if(length(unique(DATASET$Cancer)) != 1){
    VJUST_1=-1.65
    VJUST_2=-0.130
    VJUST_3=-1.25
    VJUST_4=0.275
    
    HJUST_1=0.080
    HJUST_2=0.15
    
    HJUST_3=0.080
    HJUST_4=0.15
  }
  
  if(length(unique(DATASET$Cancer)) == 1){
    VJUST_1=-1.65
    VJUST_2=-0.150
    VJUST_3=-1.25
    VJUST_4=0.275
    
    HJUST_1=-0.25
    HJUST_2=-0.6225
    HJUST_3=-0.25
    HJUST_4=-0.515
  }
  
  #plot of BAC per locus per individual
  BAC_plot <- ggplot(MAF_subset, aes(x=BAC))+ xlab('BAC') + ylab('Density') +  facet_grid(MIN_READS~TYPE) +
    
    geom_line(stat='density',size=1,alpha=1) + theme(axis.title.y = element_text(size=32)) +
    geom_area(data= AUC_coords, aes(x = x, y = y),fill='red',alpha=0.10) +
    
    geom_segment(aes(x=1/MINIMUM_READS,xend=1/MINIMUM_READS,y=0,yend=AUC_coords[which.min(abs(1/MINIMUM_READS - AUC_coords$x)),]$y),linetype='dashed',color='#4DAF4A',size=1.0) +
    
    geom_text(aes(x=1/MINIMUM_READS,y=AUC_coords[which.min(abs(1/MINIMUM_READS - AUC_coords$x)),]$y,label = paste0(1/MINIMUM_READS), vjust = -1.10,hjust=0.5),color='black',size=6.5) +
    
    #ecdf indicates proportion of samples with value below cutoff, e.g. BAC below 5; if look at MAC, can fonrim that .2% of all variants have BAC below 0.025
    #so maybe ecdf not the most accurate b/c doesn't go all the way to zero; instead, calculate actual proportion of values below these levels
    #    geom_text(aes(x=.20,y=55000,label = paste0(round(as.numeric(ecdf(MAF_subset$BAC)(1/MINIMUM_READS)-ecdf(MAF_subset$BAC)(min(AUC_coords$x))),digits=3)*100,"% of variants potentially missed\n (min. BAC = ",1/MINIMUM_READS,")"), vjust = -0.20,hjust=-0.050),size=4) + #note that x-value goes in parentheses after the ecdf fxn
    
    {if(MIN_READS == 5 | MIN_READS == 10){geom_text(aes(x=X_POS,y=Y_POS,label = paste0(round(nrow(subset(MAF_subset, BAC < 1/MINIMUM_READS))/nrow(MAF_subset),digits=3)*100,"% of variants potentially missed"), vjust = VJUST_1,hjust=HJUST_1),size=8)}} +  
    {if(MIN_READS == 5 | MIN_READS == 10){geom_text(aes(x=X_POS,y=Y_POS,label = paste0("(min. BAC = ",1/MINIMUM_READS,")"), vjust = VJUST_2,hjust=HJUST_2),size=7)}} +  
    
    {if(MIN_READS == 40){geom_text(aes(x=X_POS,y=Y_POS,label = paste0(round(nrow(subset(MAF_subset, BAC < 1/MINIMUM_READS))/nrow(MAF_subset),digits=3)*100,"% of variants potentially missed"), vjust = VJUST_3,hjust=HJUST_3),size=8)}} +  
    {if(MIN_READS ==40){geom_text(aes(x=X_POS,y=Y_POS,label = paste0("(min. BAC = ",1/MINIMUM_READS,")"), vjust = VJUST_4,hjust=HJUST_4),size=7)}} +  
    
    
    
    scale_y_continuous(expand = expansion(mult = c(.05, 0.30))) + ylab("") + 
    
    
    #scale_x_continuous(trans=pseudo_log_trans(base=2),breaks = 2^seq(1,32, by = 4)) + 
    #scale_color_manual(values=c('#00BFC4','#F8766D')) +
    scale_color_manual(values=c('red','green','blue','black','cyan')) + 
    theme(strip.text = element_text(size = 26),#,hjust=0),
          strip.text.x = element_blank(),
          axis.text.x = element_text(size=27,margin=margin(t=5,l=0,b=0,r=0)),
          axis.text.y = element_text(size=25,margin=margin(t=,l=0,b=0,r=5)),
          axis.title.y = element_blank(),
          plot.margin = unit(c(0,1,0,1), "lines"), #top, right, bottom, left
          #      axis.ticks.y=element_blank(),
          #   strip.background = element_rect(color='black', fill='white'),
          #    panel.background = element_rect(color = 'black',fill='white'),
          #    panel.grid.major = element_blank(),
          #    panel.grid.minor = element_blank(),
          #   plot.tag = element_text(size = 22,face='bold'),a
          #  plot.tag.position = c(-0.035, 0.94),
          plot.title = element_text(hjust=0.5,size=26,face='bold'),
          axis.title.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(color = 'black',fill=NA),
          strip.background = element_blank(),
          legend.position = 'right',legend.direction='vertical',
          legend.text=element_text(size=32),
          legend.title=element_blank()) +
    
    {if(MIN_READS==40){labs(title="Minimum concentration of somatic variants\n detectable at various coverage levels")}} +
    # {if(MIN_READS==40){labs(tag = "B.")}} +
    
    {if(MIN_READS!=40){labs(title="")}}
  
  return(BAC_plot)
}



MAF_protein_altering <- race_MAF_depths[grepl('Frame_Shift_Del|Frame_Shift_Ins|In_Frame_Del|In_Frame_Ins|Missense_Mutations|Nonsense_Mutation|Nonstop_Mutation|Splice_Site|Translation_Start_Site',race_MAF_depths$Variant_Classification),]
MAF_protein_altering <- subset(MAF_protein_altering, (Cancer=='BRCA' | Cancer=='LUAD' | Cancer == 'COAD' | Cancer == 'KIRC' | Cancer == 'UCEC' | Cancer == 'PRAD'))
MAF_protein_altering <- MAF_protein_altering[,c('submitter_id','race_PCA','Chromosome','Start_Position','Hugo_Symbol','t_depth','t_ref_count','t_alt_count','Variant_Classification','Cancer')]
MAF_protein_altering$BAC <- MAF_protein_altering$t_alt_count/MAF_protein_altering$t_depth

BRCA_protein_altering <- subset(MAF_protein_altering, Cancer=='BRCA')


BAC5_all <- GENERATE_BAC_PLOT(5,MAF_protein_altering,0.55,1.5)
BAC10_all <-  GENERATE_BAC_PLOT(10,MAF_protein_altering,0.55,1.5)
BAC40_all <-  GENERATE_BAC_PLOT(40,MAF_protein_altering,0.55,1.5)

BAC5_BRCA <- GENERATE_BAC_PLOT(5,BRCA_protein_altering,0.20,2)
BAC10_BRCA <-  GENERATE_BAC_PLOT(10,BRCA_protein_altering,0.20,2)
BAC40_BRCA <-  GENERATE_BAC_PLOT(40,BRCA_protein_altering,0.20,2)

BAC_fig_all <- ggarrange(BAC40_all,BAC10_all,BAC5_all,NULL,nrow=4,heights=c(1,1,1,0.05))
BAC_fig_all <- annotate_figure(BAC_fig_all,left = text_grob("Density",size=29,rot = 90),bottom = text_grob("B allele concentration (BAC)",size=29)) + theme(plot.margin = margin(0.1,0.1,0.75,0.1, "cm")) #incr marrgin at bottom

BAC_fig_BRCA <- ggarrange(BAC40_BRCA,BAC10_BRCA,BAC5_BRCA,NULL,nrow=4,heights=c(1,1,1,0.05))
BAC_fig_BRCA <- annotate_figure(BAC_fig_BRCA,left = text_grob("Density",size=29,rot = 90),bottom = text_grob("B allele concentration (BAC)",size=29)) + theme(plot.margin = margin(0.1,0.1,0.50,0.1, "cm")) #incr marrgin at bottom



######################################
## of ALT reads per locus per exome
######################################


ALT_READS_PLOT <- function(DATASET){
  
  PLOT <- ggplot(subset(DATASET,t_alt_count > 0), aes(x=t_alt_count))+ labs(title='# of reads supporting somatic\nvariant alternative alleles') + xlab('# of reads') + ylab('Number of somatic variants\n(log2 scale)') +#  facet_wrap(~race) +
    
    geom_bar(aes(x=t_alt_count),color='black',fill='cyan') +
    
    scale_y_continuous(trans=pseudo_log_trans(base=2),breaks = 2^seq(0,16, by = 2),expand=c(0.01,0.01)) +
    #scale_x_continuous(sec.axis=sec_axis(trans=~ . * (1/.)^2, name="displacement (L)"), expand = c(0, 0))
    
    #scale_x_continuous(trans=pseudo_log_trans(base=2),breaks = 2^seq(1,32, by = 4)) + 
    #scale_color_manual(values=c('#00BFC4','#F8766D')) +
    scale_color_manual(values=c('red','green','blue','black','cyan')) + 
    theme(strip.text = element_text(size = 21,face='bold'),#,hjust=0),
          axis.text.x = element_text(size=27,angle = 45,margin=margin(t=5,r=0,b=0,l=0),hjust = 1),
          axis.text.y = element_text(size=27,margin=margin(r=5)),
          axis.title.y = element_text(size=30,margin=margin(t=0,r=10,b=0,l=0)),
          axis.title.x = element_text(size=30,margin=margin(t=4,r=0,b=0,l=0)),
          #  axis.ticks= element_blank(), 
          strip.background = element_rect(color='black', fill='white'),
          panel.background = element_rect(color = 'black',fill=NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(size=28,face='bold',hjust=0.5),
          legend.position = 'right',legend.direction='vertical',
          legend.text=element_text(size=32),
          legend.title=element_blank())  + 
    scale_x_continuous(trans=pseudo_log_trans(base=2),breaks = 2^seq(0,12, by = 1),
                       expand = expansion(mult = c(.04, 0.15))) 
  
  return(PLOT)
}


variant_count_plot_all <- ALT_READS_PLOT(race_MAF_depths)
race_MAF_depths_BRCA <- subset(race_MAF_depths, Cancer=='BRCA')
variant_count_plot_BRCA <- ALT_READS_PLOT(race_MAF_depths_BRCA)




#OUTPUT FIGURE 4
pdf(paste0("/home/mayo/m187735/s212975.Wickland_Immunomics/GitHub/Racial_Disparities_in_Exome_Read_Depth/MS_plots/FIG4.pdf"),height=20,width=24)
ROW_1 <- ggarrange(MOSAIC_TUMOR,NULL,widths=c(1,0.05),font.label=list(size=32,color='black'),labels=c('A'))
ROW_3 <- ggarrange(NULL,ggarrange(variant_count_plot_BRCA,NULL,nrow=2,heights=c(1,0.3)),NULL,BAC_fig_BRCA,NULL,widths=c(0.035,0.525,0.03,0.525,0.005),nrow=1,font.label=list(size=32,color='black'),labels=c('B','','','C',''))
ggarrange(NULL,ROW_1,NULL,ROW_3,NULL,widths=c(1,0.1,1,0.05),nrow=5,heights=c(0.01,0.45,0.015,0.65,0.005))
dev.off()

#OUTPUT FIGURE 4
#pdf(paste0("/home/mayo/m187735/s212975.Wickland_Immunomics/GitHub/Racial_Disparities_in_Exome_Read_Depth/MS_plots/FIG4.pdf"),height=20,width=24)
#ggarrange(NULL,ROW_1,NULL,ROW_3,NULL,widths=c(1,0.1,1,0.05),nrow=5,heights=c(0.01,0.45,0.015,0.65,0.005))
#dev.off()

#OUTPUT FIGURE S9
#pdf(paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/GitHub/Racial_Disparities_in_Exome_Read_Depth/MS_plots/FIGS6.pdf'),height=17,width=13)
png(paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/GitHub/Racial_Disparities_in_Exome_Read_Depth/MS_plots/FIGS9.png'),height=17,width=13,res = 1400,units = 'in')
ggarrange(NULL,ggarrange(variant_count_plot_all,NULL,nrow=1,widths=c(1,0.05)),NULL,NULL,NULL,ggarrange(BAC_fig_all,NULL,nrow=1,widths=c(1,0.01)),nrow=3,ncol=2,widths=c(0.05,1),heights=c(0.75,0.05,1.50),labels=c('A.','','','','B.',''),font.label = list(size = 26,color = 'black'))
dev.off()

#####################################
#QUALITY AT SOMATIC VARIANT POSITIONS IN MAF
######################################

#Create dataframe to hold MAF data for all cancers
protected_MAF_all <- data.frame()

CANCER_LIST <- c('BRCA','UCEC','PRAD','LUAD','COAD','KIRP','KIRC')

for (CANCER_NAME in CANCER_LIST){
  print(CANCER_NAME)
  protected_MAF_cancer <- read.table(paste0('/research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/processing/TCGA/MAFs/protected_hg38/',CANCER_NAME,'_protected.maf'),header=TRUE,quote="",sep='\t')[,c('Tumor_Sample_Barcode','Chromosome','Start_Position','End_Position','Reference_Allele','Tumor_Seq_Allele1','Tumor_Seq_Allele2','n_depth','n_ref_count','n_alt_count','t_depth','t_ref_count','t_alt_count','vcf_info','vcf_format')]
  
  protected_MAF_cancer$submitter_id <- sub("^([^-]*-[^-]*-[^-]*).*", "\\1",protected_MAF_cancer$Tumor_Sample_Barcode)
  # names(protected_MAF_cancer)[2:3] <- c('tumor_tissue_depth_protected_MAF','normal_tissue_depth_protected_MAF')
  protected_MAF_cancer$Cancer <- CANCER_NAME
  protected_MAF_all <- rbind(protected_MAF_all,protected_MAF_cancer)
}
rm(protected_MAF_cancer)



#extract somatic MAF positions from protected MAF
somatic_with_vcf_and_normal_info <- merge(race_MAF_depths,protected_MAF_all,by=c('Cancer','Tumor_Sample_Barcode','Chromosome','Start_Position','End_Position','Reference_Allele','Tumor_Seq_Allele1','Tumor_Seq_Allele2','t_depth','t_ref_count','t_alt_count','submitter_id','n_depth'))

somatic_with_vcf_and_normal_info$TLOD <-sub(".*NLOD.*\\;","\\1",somatic_with_vcf_and_normal_info$vcf_info)
somatic_with_vcf_and_normal_info$TLOD <-sub("TLOD=","",somatic_with_vcf_and_normal_info$TLOD)

somatic_with_vcf_and_normal_info$NLOD <- sub(".*NLOD","\\1",somatic_with_vcf_and_normal_info$vcf_info)
somatic_with_vcf_and_normal_info$NLOD <- sub(";.*","\\1",somatic_with_vcf_and_normal_info$NLOD)
somatic_with_vcf_and_normal_info$NLOD <- sub("=","",somatic_with_vcf_and_normal_info$NLOD)


#if any barcode-allele-depth has multiple VCFs, keep the ony with highest NLOD/TLOD
#identify all duplicated rows (not just 2nd occurrence)
#in some cases multiple VCFs per sample / allele -- keep one with  largest NLOD/TLOD (affects ~110 samples)
somatic_with_vcf_and_normal_info <- somatic_with_vcf_and_normal_info %>% group_by(Cancer,Tumor_Sample_Barcode,Chromosome,Start_Position,End_Position,Reference_Allele,Tumor_Seq_Allele1,Tumor_Seq_Allele2,t_depth,t_ref_count,t_alt_count,submitter_id) %>% top_n(1,abs(as.numeric(NLOD)+as.numeric(TLOD)))

#get rid of duplicates
somatic_with_vcf_and_normal_info <- unique(somatic_with_vcf_and_normal_info)

#subset to keep only samples analyzed in previous analyses
somatic_with_vcf_and_normal_info$submitter_id_med <-  gsub("^([^-]*-[^-]*-[^-]*-[^-]*).*", "\\1",somatic_with_vcf_and_normal_info$Tumor_Sample_Barcode)

somatic_with_vcf_and_normal_info <- somatic_with_vcf_and_normal_info[somatic_with_vcf_and_normal_info$submitter_id_med %in% MAPPED_READS_FOR_PLOTTING$submitter_id_med,]

#melt the dataframe so that can have separate tumor and normal facets
somatic_with_vcf_and_normal_info_melted <- reshape2::melt(somatic_with_vcf_and_normal_info, id.vars=names(somatic_with_vcf_and_normal_info)[!grepl("NLOD|TLOD",colnames(somatic_with_vcf_and_normal_info))])

#arrange cancer order
somatic_with_vcf_and_normal_info_melted$Cancer <- factor(somatic_with_vcf_and_normal_info_melted$Cancer,levels=c('BRCA','UCEC','PRAD','LUAD','COAD','KIRP','KIRC'),ordered=TRUE)


#####################################
#QUALITY AT SOMATIC VARIANT POSITIONS IN MAF
######################################

#Create dataframe to hold MAF data for all cancers
protected_MAF_all <- data.frame()

CANCER_LIST <- c('BRCA','UCEC','PRAD','LUAD','COAD','KIRP','KIRC')

for (CANCER_NAME in CANCER_LIST){
  print(CANCER_NAME)
  protected_MAF_cancer <- read.table(paste0('/research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/processing/TCGA/MAFs/protected_hg38/',CANCER_NAME,'_protected.maf'),header=TRUE,quote="",sep='\t')[,c('Tumor_Sample_Barcode','Chromosome','Start_Position','End_Position','Reference_Allele','Tumor_Seq_Allele1','Tumor_Seq_Allele2','n_depth','n_ref_count','n_alt_count','t_depth','t_ref_count','t_alt_count','vcf_info','vcf_format')]
  
  protected_MAF_cancer$submitter_id <- sub("^([^-]*-[^-]*-[^-]*).*", "\\1",protected_MAF_cancer$Tumor_Sample_Barcode)
  # names(protected_MAF_cancer)[2:3] <- c('tumor_tissue_depth_protected_MAF','normal_tissue_depth_protected_MAF')
  protected_MAF_cancer$Cancer <- CANCER_NAME
  protected_MAF_all <- rbind(protected_MAF_all,protected_MAF_cancer)
}
rm(protected_MAF_cancer)



#extract somatic MAF positions from protected MAF
somatic_with_vcf_and_normal_info <- merge(race_MAF_depths,protected_MAF_all,by=c('Cancer','Tumor_Sample_Barcode','Chromosome','Start_Position','End_Position','Reference_Allele','Tumor_Seq_Allele1','Tumor_Seq_Allele2','t_depth','t_ref_count','t_alt_count','submitter_id','n_depth'))

somatic_with_vcf_and_normal_info$TLOD <-sub(".*NLOD.*\\;","\\1",somatic_with_vcf_and_normal_info$vcf_info)
somatic_with_vcf_and_normal_info$TLOD <-sub("TLOD=","",somatic_with_vcf_and_normal_info$TLOD)

somatic_with_vcf_and_normal_info$NLOD <- sub(".*NLOD","\\1",somatic_with_vcf_and_normal_info$vcf_info)
somatic_with_vcf_and_normal_info$NLOD <- sub(";.*","\\1",somatic_with_vcf_and_normal_info$NLOD)
somatic_with_vcf_and_normal_info$NLOD <- sub("=","",somatic_with_vcf_and_normal_info$NLOD)


#if any barcode-allele-depth has multiple VCFs, keep the ony with highest NLOD/TLOD
#identify all duplicated rows (not just 2nd occurrence)
#in some cases multiple VCFs per sample / allele -- keep one with  largest NLOD/TLOD (affects ~110 samples)
somatic_with_vcf_and_normal_info <- somatic_with_vcf_and_normal_info %>% group_by(Cancer,Tumor_Sample_Barcode,Chromosome,Start_Position,End_Position,Reference_Allele,Tumor_Seq_Allele1,Tumor_Seq_Allele2,t_depth,t_ref_count,t_alt_count,submitter_id) %>% top_n(1,abs(as.numeric(NLOD)+as.numeric(TLOD)))

#get rid of duplicates
somatic_with_vcf_and_normal_info <- unique(somatic_with_vcf_and_normal_info)

#subset to keep only samples analyzed in previous analyses
somatic_with_vcf_and_normal_info$submitter_id_med <-  gsub("^([^-]*-[^-]*-[^-]*-[^-]*).*", "\\1",somatic_with_vcf_and_normal_info$Tumor_Sample_Barcode)

somatic_with_vcf_and_normal_info <- somatic_with_vcf_and_normal_info[somatic_with_vcf_and_normal_info$submitter_id_med %in% MAPPED_READS_FOR_PLOTTING$submitter_id_med,]

#melt the dataframe so that can have separate tumor and normal facets
somatic_with_vcf_and_normal_info_melted <- reshape2::melt(somatic_with_vcf_and_normal_info, id.vars=names(somatic_with_vcf_and_normal_info)[!grepl("NLOD|TLOD",colnames(somatic_with_vcf_and_normal_info))])

#arrange cancer order
somatic_with_vcf_and_normal_info_melted$Cancer <- factor(somatic_with_vcf_and_normal_info_melted$Cancer,levels=c('BRCA','UCEC','PRAD','LUAD','COAD','KIRP','KIRC'),ordered=TRUE)

#calculte couns, avg qual
Counts_QUAL <- aggregate(value ~ variable + Cancer + race_PCA, somatic_with_vcf_and_normal_info_melted, length)
Means_QUAL <- aggregate(as.numeric(as.character(value)) ~ variable + Cancer + race_PCA, somatic_with_vcf_and_normal_info_melted, mean)
names(Means_QUAL)[4] <- 'value'



#create labels for the annotated text below each facet row
TLOD_label_counts <- data.frame(race_PCA='White',variable=c('TLOD'), read_depth=c(14000000),Cancer=c('BRCA','UCEC','PRAD','LUAD','COAD','KIRC'),label=c('n','','','','',''))
TLOD_label_counts$Cancer <- factor(TLOD_label_counts$Cancer, levels=c('BRCA','UCEC','PRAD','LUAD','COAD','KIRC'))
NLOD_label_counts <- data.frame(race_PCA='White',variable=c('NLOD'), read_depth=c(14000000),Cancer=c('BRCA','UCEC','PRAD','LUAD','COAD','KIRC'),label=c('n','','','','',''))
NLOD_label_counts$Cancer <- factor(NLOD_label_counts$Cancer, levels=c('BRCA','UCEC','PRAD','LUAD','COAD','KIRC'))

TLOD_label_mean <- data.frame(race_PCA='White',variable=c('TLOD'), read_depth=c(14000000),Cancer=c('BRCA','UCEC','PRAD','LUAD','COAD','KIRC'),label=c('bar(x)','','','','',''))
TLOD_label_mean$Cancer <- factor(TLOD_label_mean$Cancer, levels=c('BRCA','UCEC','PRAD','LUAD','COAD','KIRC'))
NLOD_label_mean <- data.frame(race_PCA='White',variable=c('NLOD'), read_depth=c(14000000),Cancer=c('BRCA','UCEC','PRAD','LUAD','COAD','KIRC'),label=c('bar(x)','','','','',''))
NLOD_label_mean$Cancer <- factor(NLOD_label_mean$Cancer, levels=c('BRCA','UCEC','PRAD','LUAD','COAD','KIRC'))


#plot
MAF_quality_by_race <- ggplot(somatic_with_vcf_and_normal_info_melted, aes(x=factor(race_PCA),y=as.numeric(value),fill=race_PCA)) + xlab("") + ylab("Log-likelihood ratio per position (log2 scale)") + facet_grid(variable~Cancer, labeller = label_wrap_gen(width=19)) + labs(title='Confidence of somatic mutation call per variant position') + 
   geom_violin() +  stat_summary(geom='crossbar',fun='median',color='white',fatten=2,width=0.35) + guides(fill = guide_legend(override.aes = list(linetype = 0))) +  
  theme(strip.text.x = element_text(size = 24),#,hjust=0),
        strip.text.y = element_text(size = 22),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=23,margin=margin(r=7)),
        axis.title.y = element_text(size=24,margin=margin(t=0,r=10,b=0,l=0)),
        axis.ticks.x = element_blank(), 
        plot.title = element_text(hjust=0.5,vjust=4,size=28,face='bold'),
        plot.tag = element_text(size = 32),
        plot.tag.position = c(-0.01, 0.98),
        plot.margin = unit(c(2,0,0,0), "lines"), #top, right, bottom, left
        strip.background = element_blank(),
        panel.background = element_rect(color = 'black',fill=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.y=unit(5,"lines"), #to control distance betweeen facet rows
        legend.position = 'bottom',legend.direction='horizontal',
        legend.text=element_text(size=32),
        legend.background = element_rect(linetype='solid', colour='black'),
        legend.title=element_blank()) + #labs(tag = "D.") +
  #geom_point() + 
  # ylim(-50000000,450000000) + 
  
  scale_fill_manual(values=c('#EF5350','#42A5F5')) + #MANUALLY SET FILL COLORS 
  
  
  scale_y_continuous(trans=pseudo_log_trans(base=2),breaks = 2^seq(1,14, by = 2),
                     expand = expansion(mult = c(0, 0.175))) +
  
  
  #means comparisons
  stat_compare_means(size=8,method='wilcox.test',mapping=aes(label=..p.signif..),
                     comparisons = list(c("Black","White")),bracket.size = NA,tip.length = 0,color='black',
                     #label.y=400000000, #adjust position of p-value
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("p<0.0001", "p<0.001", "p<0.01", "p<0.05", ""))) +
  
  coord_cartesian(ylim=c(0,NA),clip='off') + 
  
  #Add counts means
  geom_text(data = subset(Counts_QUAL, variable == "TLOD"),aes(label=formatC(value,format='f',big.mark=',',digits=0), y=-3,vjust=0.2),size=7) +
  geom_text(data = subset(Counts_QUAL, variable == "NLOD"),aes(label=formatC(value,format='f',big.mark=',',digits=0), y=-2.5,vjust=0.2),size=7) +
  
  geom_text(data = subset(Means_QUAL, variable == "TLOD"),aes(label=formatC(value,format='f',big.mark=',',digits=1), y=-9.5,vjust=0.2),size=7) +
  geom_text(data = subset(Means_QUAL, variable == "NLOD"),aes(label=formatC(value,format='f',big.mark=',',digits=1), y=-8.5,vjust=0.2),size=7) +
  
  #add text annotations
  geom_text(data= TLOD_label_counts,aes(label=label,y=-3.25,vjust=0,hjust=4.75,fontface="bold"),size=7) +
  geom_text(data= TLOD_label_mean,aes(label=label,y=-8.5,vjust=0.4,hjust=5.70,fontface="bold"),size=7,parse = TRUE) +
  geom_text(data= TLOD_label_mean,aes(label=label,y=-8.5,vjust=0.4,hjust=5.70,fontface="bold"),size=7,parse=TRUE) + #}} +
  geom_text(data= TLOD_label_mean,aes(label=label,y=-8.5,vjust=0.4,hjust=5.70,fontface="bold"),size=7.05,parse=TRUE) + #}} +
  geom_text(data= TLOD_label_mean,aes(label=label,y=-8.5,vjust=0.4,hjust=5.70,fontface="bold"),size=7.10,parse=TRUE) + #}} +
  geom_text(data= TLOD_label_mean,aes(label=label,y=-8.5,vjust=0.4,hjust=5.70,fontface="bold"),size=7.12,parse=TRUE) + #}} +
  geom_text(data= TLOD_label_mean,aes(label=label,y=-8.5,vjust=0.4,hjust=5.70,fontface="bold"),size=7.15,parse=TRUE) + #}} +
  geom_text(data= TLOD_label_mean,aes(label=label,y=-8.5,vjust=0.4,hjust=5.70,fontface="bold"),size=7.20,parse=TRUE) + #}} +
  geom_text(data= TLOD_label_mean,aes(label=label,y=-8.5,vjust=0.4,hjust=5.75,fontface="bold"),size=7.18,parse=TRUE) + #}} +
  
  
  #  {if(EXOME_OR_RNASEQ=='exome'){
  geom_text(data= NLOD_label_counts,aes(label=label,y=-2.75,vjust=0,hjust=4.75,fontface="bold"),size=7) +
  geom_text(data= NLOD_label_mean,aes(label=label,y=-8,vjust=0.2,hjust=5.70,fontface="bold"),size=7,parse=TRUE) + #}} +
  geom_text(data= NLOD_label_mean,aes(label=label,y=-8,vjust=0.2,hjust=5.70,fontface="bold"),size=7.05,parse=TRUE) + #}} +
  geom_text(data= NLOD_label_mean,aes(label=label,y=-8,vjust=0.2,hjust=5.70,fontface="bold"),size=7.10,parse=TRUE) + #}} +
  geom_text(data= NLOD_label_mean,aes(label=label,y=-8,vjust=0.2,hjust=5.70,fontface="bold"),size=7.12,parse=TRUE) + #}} +
  geom_text(data= NLOD_label_mean,aes(label=label,y=-8,vjust=0.2,hjust=5.70,fontface="bold"),size=7.15,parse=TRUE) + #}} +
  geom_text(data= NLOD_label_mean,aes(label=label,y=-8,vjust=0.2,hjust=5.70,fontface="bold"),size=7.20,parse=TRUE) + #}} +
  geom_text(data= NLOD_label_mean,aes(label=label,y=-8,vjust=0.2,hjust=5.70,fontface="bold"),size=7.18,parse=TRUE)  #}} +
  
  
  
  
  
  #OUTPUT FIGURE S8
  png(paste0("/home/mayo/m187735/s212975.Wickland_Immunomics/GitHub/Racial_Disparities_in_Exome_Read_Depth/MS_plots/FIGS8.png"),height=10,width=19,units='in',res=500)
ggarrange(MAF_quality_by_race,NULL,nrow=2,common.legend = TRUE,legend = 'bottom',align='hv',heights=c(0.8,0.055),labels=c('',''),font.label = list(size = 22, color = "black")) + theme(plot.margin = margin(b=0.5,l=0.1,t=0.1,r=0.1, "cm")) #incr marrgin at bottom
dev.off()









#formatC(length(x), format="f", big.mark = ",", digits=0)


#max and min TLOD/NLOD for each cancer
aggregate(TLOD ~  Cancer + race, somatic_with_vcf_and_normal_info, max)
aggregate(TLOD ~  Cancer + race, somatic_with_vcf_and_normal_info, min)
aggregate(NLOD ~  Cancer + race, somatic_with_vcf_and_normal_info, max)
aggregate(NLOD ~  Cancer + race, somatic_with_vcf_and_normal_info, min)


