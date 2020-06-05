
library(jsonlite,warn.conflicts = FALSE)
library(tidyr,warn.conflicts = FALSE)
library(reshape2,warn.conflicts = FALSE)
library(EnvStats,warn.conflicts = FALSE)
library(dplyr,warn.conflicts = FALSE)
library(plyr)
library(readr)
library(cowplot,warn.conflicts = FALSE)
library(factoextra,warn.conflicts = FALSE)
library(ggpubr,warn.conflicts = FALSE)
library(ggfortify,warn.conflicts = FALSE)
#library(pca3d,warn.conflicts = FALSE)
library(corrplot,warn.conflicts = FALSE)
library(ggplot2,warn.conflicts = FALSE)
library(scales,warn.conflicts = FALSE)
library(stringr)
setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
theme_set(theme_grey())


CANCER_LIST <- c('BRCA','LUAD','UCEC','KIRC','PRAD','COAD','KIRP','GBM','HNSC')
#GENE_LIST <- c('TP53','PIK3CA','MUC16','USH2A','TTN')
GENE <- 'TP53'


###############################################
#0. LOAD AND PROCESS DATA
############################################### 

combined_data_master <- data.frame()
cdriver_master <- data.frame()
all_bams_WhiteExclusive_mutations_from_MAF_master <- data.frame()

#Load data
for (CANCER_NAME in CANCER_LIST){
  #load overall depths data
  combined_data <- read.table(paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/overall_depths/',CANCER_NAME,'_overall_depths_etc.txt'),header=TRUE)
  
  #load gene depths data
  cdriver <- read.table(paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/gene_depths_pairs_from_MAF/',GENE,'/',CANCER_NAME,'_',GENE,'_depths_etc.txt'),header=TRUE)
  
  #load profiling of white-exclusive-in-MAF variants in bams
  all_bams_WhiteExclusive_mutations_from_MAF <- read.table(paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/gene_depths_WhiteExlusive_from_MAF_quered_in_all_bams/',GENE,'/',CANCER_NAME,'_',GENE,'_depths_WhiteExclusive_mutations_from_MAF_queried_in_all_bams_etc.txt'),header=TRUE)
  all_bams_WhiteExclusive_mutations_from_MAF$Cancer <- CANCER_NAME
  
  #  all_bams_all_mutations <- read.table(paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/gene_depths_pairs_missing_from_VCF/',GENE,'/',CANCER_NAME,'_',names(GENE_AND_ENSEMBL_ID),'_depths_all_mutations_all_samples_etc.txt'),header=TRUE)
  #  all_bams_all_mutations$Cancer <- CANCER_NAME
  
  
  #identify only those positions that overlap between White and Black
  cdriver_white <- unique(subset(cdriver, (race == 'White'))[,c('Chromosome','Start_Position','Tumor_Seq_Ref','Tumor_Seq_Alt','Cancer')]) 
  cdriver_black <- unique(subset(cdriver, (race == 'Black'))[,c('Chromosome','Start_Position','Tumor_Seq_Ref','Tumor_Seq_Alt','Cancer')]) 
  shared_positions_cdriver <- cdriver_black[cdriver_black$Start_Position %in% cdriver_white$Start_Position,]$Start_Position
  cdriver <- cdriver[cdriver$Start_Position %in% shared_positions_cdriver,]
  
  cdriver_white <- unique(subset(cdriver, (race == 'White'))[,c('Chromosome','Start_Position','Tumor_Seq_Ref','Tumor_Seq_Alt','Cancer')]) 
  cdriver_black <- unique(subset(cdriver, (race == 'Black'))[,c('Chromosome','Start_Position','Tumor_Seq_Ref','Tumor_Seq_Alt','Cancer')]) 
  shared_positions_cdriver <- cdriver_black[cdriver_black$Start_Position %in% cdriver_white$Start_Position,]$Start_Position
  cdriver <- cdriver[cdriver$Start_Position %in% shared_positions_cdriver,]
  
  
  
  
  
  combined_data_master <- rbind(combined_data_master,combined_data)
  cdriver_master <- rbind(cdriver_master,cdriver)
  all_bams_WhiteExclusive_mutations_from_MAF_master <- rbind(all_bams_WhiteExclusive_mutations_from_MAF_master,all_bams_WhiteExclusive_mutations_from_MAF)
  
  }

#order by race
combined_data_master <- subset(combined_data_master, (race == 'White' | race == 'Black'))
combined_data_master$race <- factor(combined_data_master$race, levels = c("White","Black"),ordered=TRUE)
cdriver_master <- subset(cdriver_master, (race == 'White' | race == 'Black'))
cdriver_master$race <- factor(cdriver_master$race, levels = c("White","Black")) 
all_bams_WhiteExclusive_mutations_from_MAF_master$race <- factor(all_bams_WhiteExclusive_mutations_from_MAF_master$race, levels = c("White","Black")) 

#all_bams_WhiteExclusive_mutations_from_MAF_absent_from_VCF_master <- subset(all_bams_WhiteExclusive_mutations_from_MAF_absent_from_VCF_master,(alt_base_count_bam != 0))

#remove excess dataframes
rm(combined_data)
rm(cdriver)
rm(all_bams_WhiteExclusive_mutations_from_MAF)




###############################################
#1. TOTAL NUMBER OF READS
############################################### 
pdf("/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/plots/total_exome_reads_per_patient.pdf",height=6,width=10.5)

n_fun <- function(x){return(data.frame(y=-12000000,label=paste(length(x),"patients")))} #function for calculating n

#EXOME
ggplot(combined_data_master, aes(x=factor(Cancer),y=read_depth_exome,fill=Cancer)) + xlab("") + ylab("Number of reads per patient") +
  geom_boxplot() +  labs(title=paste("Total # of exome reads per patient",sep='')) + 
  theme(strip.text = element_text(size = 8),
        axis.text.x = element_blank(),
        axis.ticks= element_blank(), 
        plot.title = element_text(hjust=0.5,size=18,face='bold'), 
        axis.text.y = element_text(size=12), 
        axis.title.y = element_text(size=14,face='bold',margin=margin(t=0,r=10,b=0,l=0)), 
        legend.text=element_text(size=14),
        legend.title=element_blank(),legend.position = 'bottom') + geom_point() + 
  guides(fill=guide_legend(nrow=1)) + 
  #coord_cartesian(ylim = c(((-max(combined_data_master$read_depth_exome,na.rm=TRUE)*0.25)/40), max(combined_data_master$read_depth_exome,na.rm=TRUE)*1.30)) +
  stat_summary(fun.data=n_fun,geom='text',size=4) # #calculate n 

dev.off()

###############################################
#2. TOTAL NUMBER OF EXOME READS BY RACE
############################################### 

#function to calculate n and place below each box
n_fun <- function(x){return(data.frame(y=-40000000,label=paste(length(x),"patients")))} 

#plot
pdf("/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/plots/total_exome_reads_per_patient_by_race.pdf",height=8,width=15)

ggplot(combined_data_master, aes(x=factor(race),y=read_depth_exome,fill=race)) + xlab("") + ylab("# of reads per patient") + facet_wrap(~Cancer,ncol=5) + labs(title="Total # of exome reads per patient by race") + 
  geom_boxplot() + #fill=rep(c('#F8766D','#7CAE00','#00BFC4','#C77CFF'),length(CANCER_LIST))) +   
  theme(strip.text = element_text(size = 18,face='bold'),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=18),
        axis.title.y = element_text(size=20,face='bold',margin=margin(t=0,r=10,b=0,l=0)),
        axis.ticks= element_blank(), 
        plot.title = element_text(hjust=0.5,size=22,face='bold'),
        legend.position = 'bottom',legend.direction='horizontal',
        legend.text=element_text(size=24),
        legend.title=element_blank()) + 
  geom_point() + 
  ylim(-50000000,450000000) + 
  stat_summary(fun.data=n_fun,geom='text',size=5) + #calculate n
  scale_fill_manual(values=c('#F8766D','#7CAE00','#00BFC4','#C77CFF')) + #MANUALLY SET FILL COLORS
  
  #means comparisons
  stat_compare_means(size=5,method='t.test',mapping=aes(label=..p.signif..),
                     comparisons = list(c("Black","White")),
                     #label.y=400000000, #adjust position of p-value
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("p<0.0001", "p<0.001", "p<0.01", "p<0.05", "NS"))) #+

#earlier manually set colors in geom_boxplot(fill=), so must also manually set colors in legend
#  guides(fill=guide_legend(override.aes=list(shape=15,size=5,colour=c(White="#F8766D",Black="#7CAE00",Asian="#00BFC4",Unknown="#C77CFF"))))

dev.off()

#ARRANGE MULTIPLE, SEPARATE GGPLOTS, WITH VARIABLE IN NAME~~~
#First use the gg_assign and a for-loop so that each cancer's computed separately...
#all_plots <- ggarrange(plotlist=lapply(paste0("overall_exome_",CANCER_LIST),get),common.legend=TRUE,legend='bottom',nrow=4,ncol=3)


###############################################
#3. GENE READ-DEPTH BY RACE (calculated from bam) FOR PATIENTS *WITH* MUTATIONS IN THAT GENE IN VCF
############################################### 

#function for calculating n
n_fun <- function(x){return(data.frame(y=-1.75,label=paste(length(x),"mt-patient pairs")))} 

#first, melt the dataframE
from_MAF_melted <- cdriver_master[,c("submitter_id","Chromosome","Start_Position","End_Position","Variant_Type","race","Cancer","t_ref_count","t_alt_count","tumor_ref_base_count_bam","tumor_alt_base_count_bam")]
from_MAF_melted <- melt(from_MAF_melted,id.vars=c("submitter_id","Chromosome","Start_Position","End_Position","Variant_Type","race","Cancer"),na.rm = TRUE)

#change variable names

from_MAF_melted$variable <- gsub("t_ref_count","Tumor REF (MAF)",from_MAF_melted$variable)
from_MAF_melted$variable <- gsub("t_alt_count","Tumor ALT (MAF)",from_MAF_melted$variable)

from_MAF_melted$variable <- gsub("tumor_ref_base_count_bam","Tumor REF (BAM)",from_MAF_melted$variable)
from_MAF_melted$variable <- gsub("tumor_alt_base_count_bam","Tumor ALT (BAM)",from_MAF_melted$variable)


#reorder
from_MAF_melted$variable <- factor(from_MAF_melted$variable, levels = c("Tumor REF (BAM)","Tumor ALT (BAM)","Tumor REF (MAF)","Tumor ALT (MAF)"),ordered=TRUE) 


#first the tumor samples
pdf(paste0("/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/plots/",GENE,"_depth_from_MAF_by_race_MAF.pdf"),height=20,width=15)

ggplot(subset(from_MAF_melted,(variable == 'Tumor REF (MAF)' | variable == 'Tumor ALT (MAF)' | variable == 'Tumor REF (BAM)' | variable == 'Tumor ALT (BAM)')), aes(x=factor(race),y=log(value+.001,base=2),fill=race)) + xlab("") + ylab("log2(Number of reads per patient)") + facet_grid(Cancer~variable) +
  geom_boxplot() +  labs(title=paste0("Number of exome reads per patient\n for ", GENE," mutations listed in protected MAF")) + 
  theme(strip.text = element_text(size = 18,face='bold'),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=18),
        axis.title.y = element_text(size=20,face='bold',margin=margin(t=0,r=10,b=0,l=0)),
        axis.ticks= element_blank(), 
        plot.title = element_text(hjust=0.5,size=22,face='bold'),
        legend.position = 'bottom',legend.direction='horizontal',
        legend.text=element_text(size=24),
        legend.title=element_blank()) + 
  geom_point() + 
  stat_summary(fun.data=n_fun,geom='text',size=3.35) + #calculate n
  ylim(-2.5,15) +
  stat_compare_means(size=5,method='t.test',mapping=aes(label=..p.signif..),
                     comparisons = list(c("Black","White")) ,
                     #label.y=1000, #adjust position of p-value
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("p<0.0001", "p<0.001", "p<0.01", "p<0.05", "NS")))

dev.off()





###############################################
#3. GENE READ-DEPTH BY RACE (calculated from bam) FOR PATIENTS *WITHOUT* MUTATIONS IN THAT GENE IN VCF
############################################### 

#function for calculating n
n_fun <- function(x){return(data.frame(y=-1.75,label=paste(length(x),"mt-patient pairs")))} 

#first, melt the dataframE
missing_from_MAF_melted <- all_bams_WhiteExclusive_mutations_from_MAF_master[,c("submitter_id","Chromosome","Start_Position","race","Cancer","ref_base_count_bam","alt_base_count_bam")]
missing_from_MAF_melted <- melt(missing_from_MAF_melted,id.vars=c("submitter_id","Chromosome","Start_Position","race","Cancer"),na.rm = TRUE)

#change variable names
missing_from_MAF_melted$variable <- gsub("ref_base_count_bam","Tumor REF (BAM)",missing_from_MAF_melted$variable)
missing_from_MAF_melted$variable <- gsub("alt_base_count_bam","Tumor ALT (BAM)",missing_from_MAF_melted$variable)


#reorder
missing_from_MAF_melted$variable <- factor(missing_from_MAF_melted$variable, levels = c("Tumor REF (BAM)","Tumor ALT (BAM)"),ordered=TRUE) 


#first the tumor samples
pdf(paste0("/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/plots/",GENE,"_depth_missing_from_MAF_by_race_MAF.pdf"),height=20,width=15)

ggplot(subset(missing_from_MAF_melted,(variable == 'Tumor REF (BAM)' | variable == 'Tumor ALT (BAM)')), aes(x=factor(race),y=log(value+.001,base=2),fill=race)) + xlab("") + ylab("log2(Number of reads per patient)") + facet_grid(Cancer~variable) +
  geom_boxplot() +  labs(title=paste0("Number of exome reads per patient\n for White-exclusive ", GENE," mutations listed in protected MAF")) + 
  theme(strip.text = element_text(size = 18,face='bold'),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=18),
        axis.title.y = element_text(size=20,face='bold',margin=margin(t=0,r=10,b=0,l=0)),
        axis.ticks= element_blank(), 
        plot.title = element_text(hjust=0.5,size=22,face='bold'),
        legend.position = 'bottom',legend.direction='horizontal',
        legend.text=element_text(size=24),
        legend.title=element_blank()) + 
  geom_point() + 
  stat_summary(fun.data=n_fun,geom='text',size=3.35) + #calculate n
  ylim(-2.5,15) +
  stat_compare_means(size=5,method='t.test',mapping=aes(label=..p.signif..),
                     comparisons = list(c("Black","White")) ,
                     #label.y=1000, #adjust position of p-value
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("p<0.0001", "p<0.001", "p<0.01", "p<0.05", "NS")))

dev.off()


n_fun <- function(x){return(data.frame(y=-25,label=paste(length(x),"mt-patient pairs","\nmean:",round(mean(x),2),"reads","\nmedian:",round(median(x),2),"reads")))} #function for calculating n

ggplot(subset(missing_from_MAF_melted,(race=='Black' &  variable == 'Tumor ALT (BAM)')), aes(x=factor(Cancer),y=value,fill=Cancer)) + xlab("") + ylab("Number of reads per patient") +
  geom_boxplot() +  labs(title=paste0("Number of ALT reads per Black patient\n for White-exclusive ", GENE," mutations listed in protected MAF")) + 
  theme(strip.text = element_text(size = 18,face='bold'),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=18),
        axis.title.y = element_text(size=20,face='bold',margin=margin(t=0,r=10,b=0,l=0)),
        axis.ticks= element_blank(), 
        plot.title = element_text(hjust=0.5,size=22,face='bold'),
        legend.position = 'bottom',legend.direction='horizontal',
        legend.text=element_text(size=14),
        legend.title=element_blank()) + 
  guides(fill=guide_legend(nrow=1)) + 
  geom_point() + 
  stat_summary(fun.data=n_fun,geom='text',size=3.35) + #calculate n 
  ylim(-30,200) 



n_fun <- function(x){return(data.frame(y=-1,label=paste(length(x),"mt-patient pairs","\nmean:",round(mean(x),2),"reads","\nmedian:",round(median(x),2),"reads")))} #function for calculating n

ggplot(subset(missing_from_MAF_melted,(race=='Black' &  variable == 'Tumor ALT (BAM)')), aes(x=factor(Cancer),y=log(value,base=2),fill=Cancer)) + xlab("") + ylab("log2 Number of reads per patient") +
  geom_boxplot() +  labs(title=paste0("Number of ALT reads per Black patient\n for White-exclusive ", GENE," mutations listed in protected MAF")) + 
  theme(strip.text = element_text(size = 18,face='bold'),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=18),
        axis.title.y = element_text(size=20,face='bold',margin=margin(t=0,r=10,b=0,l=0)),
        axis.ticks= element_blank(), 
        plot.title = element_text(hjust=0.5,size=22,face='bold'),
        legend.position = 'bottom',legend.direction='horizontal',
        legend.text=element_text(size=14),
        legend.title=element_blank()) + 
  guides(fill=guide_legend(nrow=1)) + 
  geom_point() + 
  #stat_summary(fun.data=n_fun,geom='text',size=3.35) + #calculate n 
  ylim(0,10) 







