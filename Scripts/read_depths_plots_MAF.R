
library(jsonlite,warn.conflicts = FALSE)
library(tidyr,warn.conflicts = FALSE)
library(reshape2,warn.conflicts = FALSE)
library(EnvStats,warn.conflicts = FALSE)
library(dplyr,warn.conflicts = FALSE)
library(plyr)
library(grid)
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


CANCER_LIST <- c('BRCA')#,'LUAD','UCEC','OV','KIRC','PRAD','COAD','KIRP','GBM','HNSC')
#GENE_LIST <- c('TP53','PIK3CA','MUC16','USH2A','TTN')
GENE <- 'TP53'


###############################################
#0. LOAD AND PROCESS DATA
############################################### 

combined_data_master <- data.frame()
all_MAF_mutations_master <- data.frame()
shared_MAF_mutations_master <- data.frame()
WhiteExclusive_mutations_in_bams_VCFs_master <- data.frame()

#Load data
for (CANCER_NAME in CANCER_LIST){
  #load overall depths data
  combined_data <- read.table(paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/overall_depths/',CANCER_NAME,'_overall_depths_etc.txt'),header=TRUE)
  
  #load gene depths data
  all_MAF_mutations <- read.table(paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/gene_depths_pairs_from_MAF/',GENE,'/',CANCER_NAME,'_',GENE,'_depths_etc.txt'),header=TRUE)
  
  #load profiling of white-exclusive-in-MAF variants in bams
  WhiteExclusive_mutations_in_bams_VCFs <- read.table(paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/gene_depths_WhiteExlusive_from_MAF_quered_in_all_bams/',GENE,'/',CANCER_NAME,'_',GENE,'_depths_WhiteExclusive_mutations_from_MAF_queried_in_all_bams_and_VCFs.txt'),header=TRUE)
  WhiteExclusive_mutations_in_bams_VCFs <- WhiteExclusive_mutations_in_bams_VCFs[ , !names(WhiteExclusive_mutations_in_bams_VCFs) %in% c("chromosome")] 
  WhiteExclusive_mutations_in_bams_VCFs$Cancer <- CANCER_NAME

  
  
  #identify only those positions that overlap between White and Black
  all_MAF_mutations_white <- unique(subset(all_MAF_mutations, (race == 'White'))[,c('Chromosome','Start_Position','Tumor_Seq_Ref','Tumor_Seq_Alt','Cancer')]) 
  all_MAF_mutations_black <- unique(subset(all_MAF_mutations, (race == 'Black'))[,c('Chromosome','Start_Position','Tumor_Seq_Ref','Tumor_Seq_Alt','Cancer')]) 
  
  white_shared_with_black <- unique(all_MAF_mutations_white[all_MAF_mutations_white$Start_Position %in% all_MAF_mutations_black$Start_Position,]$Start_Position)
  black_shared_with_white <- unique(all_MAF_mutations_black[all_MAF_mutations_black$Start_Position %in% all_MAF_mutations_white$Start_Position,]$Start_Position)
  
  black_white_shared <- unique(c(white_shared_with_black, black_shared_with_white))
  shared_MAF_mutations <- all_MAF_mutations[all_MAF_mutations$Start_Position %in% black_white_shared,]
  #remove 'mutation count per patient' - no longer accurate b/c this is subset of original dataframe
  shared_MAF_mutations <- shared_MAF_mutations[-c(18)]

  
  #combine with other cancer types
  combined_data_master <- rbind(combined_data_master,combined_data)
  all_MAF_mutations_master <- rbind(all_MAF_mutations_master,all_MAF_mutations)
  shared_MAF_mutations_master <- rbind(shared_MAF_mutations_master,shared_MAF_mutations)
  WhiteExclusive_mutations_in_bams_VCFs_master <- rbind(WhiteExclusive_mutations_in_bams_VCFs_master,WhiteExclusive_mutations_in_bams_VCFs)
  
  }

#order by race
combined_data_master <- subset(combined_data_master, (race == 'White' | race == 'Black'))
combined_data_master$race <- factor(combined_data_master$race, levels = c("White","Black"),ordered=TRUE)
combined_data_master$Cancer <- factor(combined_data_master$Cancer,levels=c('BRCA','UCEC','KIRP','PRAD','COAD','LUAD','KIRC','HNSC','GBM','OV'))
all_MAF_mutations_master <- subset(all_MAF_mutations_master, (race == 'White' | race == 'Black'))
all_MAF_mutations_master$race <- factor(all_MAF_mutations_master$race, levels = c("White","Black")) 
all_MAF_mutations_master$Cancer <- factor(all_MAF_mutations_master$Cancer,levels=c('BRCA','UCEC','KIRP','PRAD','COAD','LUAD','KIRC','HNSC','GBM','OV'))
shared_MAF_mutations_master <- subset(shared_MAF_mutations_master, (race == 'White' | race == 'Black'))
shared_MAF_mutations_master$race <- factor(shared_MAF_mutations_master$race, levels = c("White","Black")) 
shared_MAF_mutations_master$Cancer <- factor(shared_MAF_mutations_master$Cancer,levels=c('BRCA','UCEC','KIRP','PRAD','COAD','LUAD','KIRC','HNSC','GBM','OV'))
WhiteExclusive_mutations_in_bams_VCFs_master$race <- factor(WhiteExclusive_mutations_in_bams_VCFs_master$race, levels = c("White","Black")) 
WhiteExclusive_mutations_in_bams_VCFs_master$Cancer <- factor(WhiteExclusive_mutations_in_bams_VCFs_master$Cancer,levels=c('BRCA','UCEC','KIRP','PRAD','COAD','LUAD','KIRC','HNSC','GBM','OV'))

#all_bams_WhiteExclusive_mutations_from_MAF_absent_from_VCF_master <- subset(all_bams_WhiteExclusive_mutations_from_MAF_absent_from_VCF_master,(alt_base_count_bam != 0))

#remove excess dataframes
rm(combined_data)
rm(all_MAF_mutations)
rm(shared_MAF_mutations)
rm(all_MAF_mutations_black, all_MAF_mutations_white)
rm(WhiteExclusive_mutations_in_bams_VCFs)



###############################################
#1. TOTAL NUMBER OF READS
############################################### 
pdf("/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/plots/total_exome_reads_per_patient.pdf",height=6,width=10.5)

n_fun <- function(x){return(data.frame(y=-12000000,label=paste(length(x),"patients")))} #function for calculating n

#EXOME
ggplot(combined_data_master, aes(x=factor(Cancer),y=read_depth_exome,fill=Cancer)) + xlab("") + ylab("# of reads per patient") +
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
  stat_compare_means(size=5.5,method='t.test',mapping=aes(label=..p.signif..),
                     comparisons = list(c("Black","White")),bracket.size = NA,tip.length = 0,color='blue',
                     #label.y=400000000, #adjust position of p-value
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("p<0.0001", "p<0.001", "p<0.01", "p<0.05", ""))) #+

#earlier manually set colors in geom_boxplot(fill=), so must also manually set colors in legend
#  guides(fill=guide_legend(override.aes=list(shape=15,size=5,colour=c(White="#F8766D",Black="#7CAE00",Asian="#00BFC4",Unknown="#C77CFF"))))

dev.off()

#ARRANGE MULTIPLE, SEPARATE GGPLOTS, WITH VARIABLE IN NAME~~~
#First use the gg_assign and a for-loop so that each cancer's computed separately...
#all_plots <- ggarrange(plotlist=lapply(paste0("overall_exome_",CANCER_LIST),get),common.legend=TRUE,legend='bottom',nrow=4,ncol=3)


###############################################
#3. GENE READ-DEPTH BY RACE FOR ALL MAF MUTATIONS
############################################### 

#function for calculating n
n_fun <- function(x){return(data.frame(y=-1.75,label=paste(length(x),"mt-patient pairs")))} 

#first, melt the dataframE
from_MAF_melted <- all_MAF_mutations_master[,c("submitter_id","Chromosome","Start_Position","End_Position","Variant_Type","race","Cancer","t_ref_count","t_alt_count","tumor_ref_base_count_bam","tumor_alt_base_count_bam")]
from_MAF_melted <- melt(from_MAF_melted,id.vars=c("submitter_id","Chromosome","Start_Position","End_Position","Variant_Type","race","Cancer"),na.rm = TRUE)

#change variable names
from_MAF_melted$variable <- gsub("t_ref_count","Tumor REF (MAF)",from_MAF_melted$variable)
from_MAF_melted$variable <- gsub("t_alt_count","Tumor ALT (MAF)",from_MAF_melted$variable)

from_MAF_melted$variable <- gsub("tumor_ref_base_count_bam","Tumor REF (BAM)",from_MAF_melted$variable)
from_MAF_melted$variable <- gsub("tumor_alt_base_count_bam","Tumor ALT (BAM)",from_MAF_melted$variable)

#reorder
from_MAF_melted$variable <- factor(from_MAF_melted$variable, levels = c("Tumor REF (MAF)","Tumor ALT (MAF)","Tumor REF (BAM)","Tumor ALT (BAM)"),ordered=TRUE) 


#plot
pdf(paste0("/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/plots/",GENE,"_depth_from_MAF_by_race.pdf"),height=20,width=15)

ggplot(subset(from_MAF_melted,(variable == 'Tumor REF (MAF)' | variable == 'Tumor ALT (MAF)' | variable == 'Tumor REF (BAM)' | variable == 'Tumor ALT (BAM)')), aes(x=factor(race),y=log(value+.001,base=2),fill=race)) + xlab("") + ylab("log2 # of reads per patient") + facet_grid(Cancer~variable) +
  geom_boxplot() +  labs(title=paste0("# of MAF & BAM exome reads per patient\n for all ", GENE," mutations reported in protected MAF")) + 
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
  stat_compare_means(size=5.5,method='t.test',mapping=aes(label=..p.signif..),
                     comparisons = list(c("Black","White")),bracket.size = NA,tip.length = 0,color='blue',
                     #label.y=400000000, #adjust position of p-value
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("p<0.0001", "p<0.001", "p<0.01", "p<0.05", ""))) #+

dev.off()



###############################################
#4. GENE READ-DEPTH BY RACE FOR MAF MUTATIONS SHARED BY BLACK AND WHITE SAMPLES
############################################### 

#function for calculating n
n_fun <- function(x){return(data.frame(y=-1.75,label=paste(length(x),"mt-patient pairs")))} 

#add column for number of mutations at locus for each race
racial_count_per_position <- shared_MAF_mutations_master %>% group_by(race,Start_Position) %>% tally(name='Racial_count_per_position')


#first, melt the dataframE
from_MAF_melted <- shared_MAF_mutations_master[,c("submitter_id","Chromosome","Start_Position","End_Position","Variant_Type","race","Cancer","t_ref_count","t_alt_count","tumor_ref_base_count_bam","tumor_alt_base_count_bam")]
from_MAF_melted <- melt(from_MAF_melted,id.vars=c("submitter_id","Chromosome","Start_Position","End_Position","Variant_Type","race","Cancer"),na.rm = TRUE)

#change variable names
from_MAF_melted$variable <- gsub("t_ref_count","Tumor REF (MAF)",from_MAF_melted$variable)
from_MAF_melted$variable <- gsub("t_alt_count","Tumor ALT (MAF)",from_MAF_melted$variable)

from_MAF_melted$variable <- gsub("tumor_ref_base_count_bam","Tumor REF (BAM)",from_MAF_melted$variable)
from_MAF_melted$variable <- gsub("tumor_alt_base_count_bam","Tumor ALT (BAM)",from_MAF_melted$variable)

#reorder
from_MAF_melted$variable <- factor(from_MAF_melted$variable, levels = c("Tumor REF (MAF)","Tumor ALT (MAF)","Tumor REF (BAM)","Tumor ALT (BAM)"),ordered=TRUE) 


#plot
pdf(paste0("/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/plots/",GENE,"_depth_from_MAF_by_race_common.pdf"),height=20,width=15)

#note that the depth in BAMs comes from the MAFs, which have cols for that info
ggplot(subset(from_MAF_melted,(variable == 'Tumor REF (MAF)' | variable == 'Tumor ALT (MAF)' | variable == 'Tumor REF (BAM)' | variable == 'Tumor ALT (BAM)')), aes(x=factor(race),y=log(value+.001,base=2),fill=race)) + xlab("") + ylab("log2 # of reads per patient") + facet_grid(Cancer~variable) +
  geom_boxplot() +  labs(title=paste0("# of MAF & BAM exome reads per patient\n for common/shared Black & White ", GENE," mutations reported in protected MAF")) + 
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
  stat_compare_means(size=5.5,method='t.test',mapping=aes(label=..p.signif..),
                     comparisons = list(c("Black","White")),bracket.size = NA,tip.length = 0,color='blue',
                     #label.y=400000000, #adjust position of p-value
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("p<0.0001", "p<0.001", "p<0.01", "p<0.05", ""))) #+

dev.off()

#For that same set of SNPs, plot by position
pdf(paste0("/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/plots/",GENE,"_depth_from_MAF_by_race_common_by_position.pdf"),height=8,width=13)

ggplot(from_MAF_melted, aes(x=as.character(Start_Position), y=as.numeric(value),color=race,label=Start_Position)) + facet_wrap(~variable,nrow=2, ncol = 2) + #treat x axis as factor so taht only values with numbers are included
  geom_smooth(aes(group=race,color=race),method='loess',se=FALSE,size=0.65) +
  geom_point(size=0.85) +
  
  theme(strip.text = element_text(size = 18,face='bold'),
        axis.text.x = element_text(vjust=0.5,size=12,face='bold',angle=90),
        axis.text.y = element_text(size=18),
        axis.title.y = element_text(size=20,face='bold',margin=margin(t=0,r=10,b=0,l=0)),
        axis.title.x = element_text(size=20,face='bold',margin=margin(t=10,r=0,b=0,l=0)),
        axis.ticks= element_blank(), 
        plot.title = element_text(hjust=0.5,size=22,face='bold'),
        legend.position = 'bottom',legend.direction='horizontal',
        legend.text=element_text(size=24),
        legend.background = element_rect(linetype='solid', colour='black'),
        legend.title=element_blank()) + 
  ggtitle("Read depth along chromosomal position in TP53") +
  xlab("Position") + ylab("# of reads per patient") 
dev.off()

###############################################
#5. GENE READ-DEPTH BY RACE FOR WHITE-EXCLUSIVE MAF MUTATIONS: LOOK UP LOCI IN BAM
############################################### 

#function for calculating n
n_fun <- function(x){return(data.frame(y=-1.75,label=paste(length(x),"mt-patient pairs")))} 

#first, melt the dataframE
from_MAF_melted <- WhiteExclusive_mutations_in_bams_VCFs_master[,c("submitter_id","Chromosome","Start_Position","End_Position","Variant_Type","race","Cancer","alt_allele_concent_bam","ref_base_count_bam","alt_base_count_bam","AD_ref_VCF","AD_alt_VCF")]

from_MAF_melted <- melt(from_MAF_melted,id.vars=c("submitter_id","Chromosome","Start_Position","End_Position","Variant_Type","race","Cancer","alt_allele_concent_bam"),na.rm = TRUE)


#change variable names
from_MAF_melted$variable <- gsub("AD_ref_VCF","Tumor REF (VCF)",from_MAF_melted$variable)
from_MAF_melted$variable <- gsub("AD_alt_VCF","Tumor ALT (VCF)",from_MAF_melted$variable)

from_MAF_melted$variable <- gsub("ref_base_count_bam","Tumor REF (BAM)",from_MAF_melted$variable)
from_MAF_melted$variable <- gsub("alt_base_count_bam","Tumor ALT (BAM)",from_MAF_melted$variable)

#reorder
from_MAF_melted$variable <- factor(from_MAF_melted$variable, levels = c("Tumor REF (BAM)","Tumor ALT (BAM)","Tumor REF (VCF)","Tumor ALT (VCF)"),ordered=TRUE) 


#plot
pdf(paste0("/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/plots/",GENE,"_depth_from_MAF_by_race_WhiteExclusive.pdf"),height=20,width=15)

ggplot(subset(from_MAF_melted,(variable == 'Tumor REF (BAM)' | variable == 'Tumor ALT (BAM)')), aes(x=factor(race),y=log(value+.001,base=2),fill=race)) + xlab("") + ylab("log2 # of reads per patient)") + facet_grid(Cancer~variable) +
  geom_boxplot() +  labs(title=paste0("# of BAM exome reads per patient\n for White-exclusive ", GENE," mutations reported in protected MAF")) + 
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
  stat_compare_means(size=5.5,method='t.test',mapping=aes(label=..p.signif..),
                     comparisons = list(c("Black","White")),bracket.size = NA,tip.length = 0,color='blue',
                     #label.y=400000000, #adjust position of p-value
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("p<0.0001", "p<0.001", "p<0.01", "p<0.05", ""))) #+

dev.off()


###################################
#6. ALT ALLELE CONCENTRATION
###################################

#function for calculating n
n_fun <- function(x){return(data.frame(y=-0.15,label=paste(length(x),"mt-patient pairs")))} 


#plot
pdf(paste0("/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/plots/",GENE,"_ALT_allele_concent_for_WhiteExclusive.pdf"),height=18,width=10)

ggplot(subset(from_MAF_melted,(variable == 'Tumor REF (BAM)' | variable == 'Tumor ALT (BAM)')), aes(x=factor(race),y=alt_allele_concent_bam,fill=race)) + xlab("") + ylab("ALT allele concentration") + facet_grid(Cancer~variable) +
  geom_violin() +  labs(title=paste0("ALT allele concentration in BAMs \nfor White-exclusive ", GENE," mutations reported in protected MAF")) + 
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
  ylim(-0.25,1.25) + 
  stat_summary(fun.data=n_fun,geom='text',size=3.35) + #calculate n
  stat_compare_means(size=5.5,method='t.test',mapping=aes(label=..p.signif..),
                     comparisons = list(c("Black","White")),bracket.size = NA,tip.length = 0,color='blue',
                     #label.y=400000000, #adjust position of p-value
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("p<0.0001", "p<0.001", "p<0.01", "p<0.05", ""))) #+
dev.off()





###################################
#7. NUMBER OF ALT READS PER BLACK PATIENT AT WHITE-EXCLUSIVE MAF LOCI: RAW COUNTS
###################################
#function for calculating n
n_fun <- function(x){return(data.frame(y=-0.125,label=paste(length(x),"mt-patient pairs")))} 

#plot
pdf(paste0("/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/plots/",GENE,"_ALT_depth_in_black__for_WhiteExclusive_cts.pdf"),height=10,width=20)

n_fun <- function(x){return(data.frame(y=-17,label=paste(length(x),"mt-patient pairs","\nmean:",round(mean(x),2),"reads","\nmedian:",round(median(x),2),"reads")))} #function for calculating n
ggplot(subset(from_MAF_melted,(race=='Black' &  (variable == 'Tumor ALT (BAM)'))), aes(x=factor(Cancer),y=value,fill=Cancer)) + xlab("") + ylab("# of reads per patient") +
  geom_boxplot() +  labs(title=paste0("# of ALT BAM reads per Black patient\n for White-exclusive ", GENE," mutations listed in protected MAF")) + 
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
  ylim(-20,200) 

dev.off()


###################################
#8. NUMBER OF ALT READS PER BLACK PATIENT AT WHITE-EXCLUSIVE MAF LOCI: LOG TRANSFORMED
###################################
#function for calculating n
n_fun <- function(x){return(data.frame(y=-1,label=paste(length(x),"mt-patient pairs")))} 
#n_fun <- function(x){return(data.frame(y=-17,label=paste(length(x),"mt-patient pairs","\nmean:",round(mean(x),2),"reads","\nmedian:",round(median(x),2),"reads")))} #function for calculating n

#plot
pdf(paste0("/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/plots/",GENE,"_ALT_depth_in_black__for_WhiteExclusive_log2.pdf"),height=10,width=20)


ggplot(subset(from_MAF_melted,(race=='Black' &  (variable == 'Tumor ALT (BAM)'))), aes(x=factor(Cancer),y=log(value+.001,base=2),fill=Cancer)) + xlab("") + ylab("Log2 # of reads per patient") +
  geom_boxplot() +  labs(title=paste0("# of ALT reads per Black patient\n for White-exclusive ", GENE," mutations listed in protected MAF")) + 
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
  stat_summary(fun.data=n_fun,geom='text',size=3.35) #+ #calculate n 
  ylim(-12,12) 

dev.off()



###############################################
#8. TOTAL NUMBER OF READS PER PATIENT BY SOURCE SITE AND RACE
############################################### 

#function for calculating n
n_fun <- function(x){return(data.frame(y=0,label=paste0("n=",length(x))))} 

#subset by cancer
source_site <- combined_data_master

#change variable names
source_site$Source_Site_exome <- gsub("_"," ",source_site$Source_Site_exome)

#calcualte number of black, white samples per site
source_site_counts <- source_site[,c('submitter_id','race','Source_Site_exome','read_depth_exome','Cancer')] %>%
  #  tidyr::complete(sample_name, HLA_type,peptide_class, fill = list(number_of_unique_alleles = 0)) %>% #retain zero-counts
  dplyr::group_by(Cancer,Source_Site_exome,race,.drop=FALSE) %>%
  dplyr::summarise(number_of_samples=n()) 

#'spread' the data and only keep centers with at least 5 black and at least 5 white samples
source_site_counts <- source_site_counts %>% spread(race, number_of_samples)
source_site_counts <- subset(source_site_counts, (Black >= 5 & White >= 5))

source_site <- source_site[paste(source_site$Cancer, source_site$Source_Site_exome) %in% paste(source_site_counts$Cancer, source_site_counts$Source_Site_exome),]

#  pdf(paste0("/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/plots/total_exome_reads_by_source_site_",CANCER_NAME,".pdf"),height=9,width=11)

for (CANCER_NAME in CANCER_LIST){  
  
  source_site_cancer <- subset(source_site, (Cancer==CANCER_NAME))
  source_site_cancer <- droplevels(source_site_cancer)
  
  if(length(unique(source_site_cancer$Source_Site_exome))==9){SCALE=0}
  if(length(unique(source_site_cancer$Source_Site_exome))==7){SCALE=7.8}
  if(length(unique(source_site_cancer$Source_Site_exome))==4){SCALE=19.85}
  if(length(unique(source_site_cancer$Source_Site_exome))==3){SCALE=23.75}
  if(length(unique(source_site_cancer$Source_Site_exome))==2){SCALE=27.70}
  
  # source_site_cancer$Source_Site_exome <- factor(source_site_cancer$Source_Site_exome, levels=c(unique(source_site_cancer$Source_Site_exome),sprintf("Blank%s",seq(1:(10-length(unique(source_site_cancer$Source_Site_exome)))))))
  
  #  pdf('/home/mayo/m187735/by_center.pdf',height=14,width=18)
  source_site_cancer_plot <- ggplot(source_site_cancer, aes(x=factor(race), y=read_depth_exome,fill=race)) + xlab("") + ylab("Exome read depth") + facet_grid(Cancer~Source_Site_exome,scales='free',drop=FALSE,labeller = label_wrap_gen(width = 20,multi_line = TRUE)) +  
    
    geom_boxplot(outlier.shape=NA) +labs(title = "") +
    theme(plot.margin=unit(c(-0.5,SCALE,-0.5,0), "cm"),
          plot.title = element_text(hjust=0.5),
          strip.text = element_text(face="bold", size=10),
          strip.text.y= element_text(size=20,color='black',face='bold'),
          axis.title.y = element_text(size=14,margin = margin(t = 0, r = 5, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 15, r = 10, b = 5, l = 10)),
          axis.text.x = element_blank(),axis.ticks= element_blank(),
          axis.text.y = element_text(size=12),
          legend.text=element_text(size=24), legend.title=element_blank(),
          legend.background = element_rect(linetype='solid', colour='black'),
          legend.position = 'none') + 
    
    #legend.position = c(0.5,-0.04),legend.direction='horizontal',plot.margin = margin(t=5,r=15,b=15,l=5))+
    geom_point() + stat_summary(fun.data=n_fun,geom='text',size=4) + #calculate n
    
    
    stat_compare_means(size=4.5,method='t.test',mapping=aes(label=..p.signif..),
                       comparisons = list(c("Black","White")),bracket.size = NA,tip.length = 0,color='blue',
                       #label.y=400000000, #adjust position of p-value
                       symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("p<0.0001", "p<0.001", "p<0.01", "p<0.05", ""))) +
    scale_y_continuous(labels = scales::label_scientific(digits=1)) +  scale_fill_manual(values=c('#F8766D','#00BFC4')) +
    coord_cartesian(ylim=c(0, max(source_site_cancer$read_depth_exome*1.15)))
  
  
  #now convert to grob, identify which cols contain facets, and set facet width
  # source_site_cancer <- ggplot_gtable(ggplot_build(source_site_cancer))
  #  facet.cols <- source_site_cancer$layout$l[grepl("panel", source_site_cancer$layout$name)]
  # source_site_cancer$widths[facet.cols] <- unit(c(5), "cm")
  # source_site_cancer$heights[facet.cols] <- unit(c(2),"null")
  assign(paste0("plot_",CANCER_NAME), source_site_cancer_plot)
}

#combine plots and output
pdf('/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/plots/total_exome_reads_by_site.pdf',height=20,width=15)
figure <- ggarrange(plot_BRCA,plot_UCEC,plot_KIRP,plot_PRAD,plot_COAD,plot_LUAD,plot_KIRC,plot_HNSC,plot_GBM,nrow=9,common.legend = TRUE, legend="bottom")
annotate_figure(figure,top = text_grob("Total # of exome reads per sampling site for sites contributing at least 5 Black and 5 White samples",face = "bold", size = 20))

dev.off()















cdriver <- add_count(cdriver, submitter_id)








