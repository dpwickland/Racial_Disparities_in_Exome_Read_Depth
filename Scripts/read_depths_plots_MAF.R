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



GENE='TP53'

CANCER_LIST <- c('BRCA','LUAD','UCEC','COAD','PRAD','KIRC')

cdriver_master <- data.frame()
#Load data
for (CANCER_NAME in CANCER_LIST){

  #load gene depths data
  cdriver <- read.table(paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/gene_depths_pairs_from_MAF/',GENE,'/',CANCER_NAME,'_',GENE,'_depths_etc.txt'),header=TRUE)
  
  cdriver_master <- rbind(cdriver_master,cdriver)

}

#order by race
cdriver_master$race <- factor(cdriver_master$race, levels = c("White","Black","Asian","Unknown")) 

#remove excess dataframes
rm(cdriver)

#subset white and black samples
cdriver_master <- subset(cdriver_master, (race == 'White' | race == 'Black'))

###############################################
#3. GENE READ-DEPTH BY RACE (calculated from bam) FOR PATIENTS *WITH* MUTATIONS IN THAT GENE IN VCF
############################################### 



#function for calculating n
n_fun <- function(x){return(data.frame(y=-1.75,label=paste(length(x),"mutation-\npatient pairs")))} 

#first, melt the dataframE
from_VCF_melted <- cdriver_master[,c("submitter_id","Chromosome","Start_Position","End_Position","Variant_Type","race","Cancer","t_depth","t_ref_count","t_alt_count","n_depth","n_ref_count","n_alt_count","variant_read_depth_bam_MPILEUP","tumor_ref_base_count_bam","tumor_alt_base_count_bam")]
from_VCF_melted <- melt(from_VCF_melted,id.vars=c("submitter_id","Chromosome","Start_Position","End_Position","Variant_Type","race","Cancer"),na.rm = TRUE)

#change variable names
from_VCF_melted$variable <- gsub("t_depth","Tumor REF+ALT (MAF)",from_VCF_melted$variable)
from_VCF_melted$variable <- gsub("t_ref_count","Tumor REF (MAF)",from_VCF_melted$variable)
from_VCF_melted$variable <- gsub("t_alt_count","Tumor ALT (MAF)",from_VCF_melted$variable)

from_VCF_melted$variable <- gsub("n_depth","Normal REF+ALT (MAF)",from_VCF_melted$variable)
from_VCF_melted$variable <- gsub("n_ref_count","Normal REF (MAF)",from_VCF_melted$variable)
from_VCF_melted$variable <- gsub("n_alt_count","Normal ALT (MAF)",from_VCF_melted$variable)

from_VCF_melted$variable <- gsub("variant_read_depth_bam_MPILEUP","Tumor REF+ALT (BAM)",from_VCF_melted$variable)
from_VCF_melted$variable <- gsub("tumor_ref_base_count_bam","Tumor REF (BAM)",from_VCF_melted$variable)
from_VCF_melted$variable <- gsub("tumor_alt_base_count_bam","Tumor ALT (BAM)",from_VCF_melted$variable)


#reorder
from_VCF_melted$variable <- factor(from_VCF_melted$variable, levels = c("Tumor REF+ALT (BAM)","Tumor REF (BAM)","Tumor ALT (BAM)","Tumor REF+ALT (MAF)","Tumor REF (MAF)","Tumor ALT (MAF)","Normal REF+ALT (MAF)","Normal REF (MAF)","Normal ALT (MAF)"),ordered=TRUE) 


#first the tumor samples
pdf(paste0("/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/plots/",GENE,"_depth_from_MAF_by_race_MAF.pdf"),height=14,width=13)

ggplot(subset(from_VCF_melted,(variable == 'Tumor REF+ALT (MAF)' | variable == 'Tumor REF (MAF)' | variable == 'Tumor ALT (MAF)')), aes(x=factor(race),y=log(value+.001,base=2),fill=race)) + xlab("") + ylab("log2(Number of reads per patient)") + facet_grid(Cancer~variable) +
  geom_boxplot() +  labs(title=paste0("Exome read depth of mutations in ", GENE,"\nfor patients WITH these ",GENE," mutations in MAF\n(mutations identified from MAF; read depth reported in MAF)")) + 
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
  stat_summary(fun.data=n_fun,geom='text',size=3.5) + #calculate n
  ylim(-2.5,13) +
  stat_compare_means(size=5,method='t.test',mapping=aes(label=..p.signif..),
                     comparisons = list(c("Black","White")) ,
                     #label.y=1000, #adjust position of p-value
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("p<0.0001", "p<0.001", "p<0.01", "p<0.05", "NS")))
dev.off()


pdf(paste0("/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/plots/",GENE,"_depth_from_MAF_by_race_BAM.pdf"),height=14,width=13)

#for now, do SNPs only
ggplot(subset(from_VCF_melted,(Variant_Type == 'SNP' & (variable == 'Tumor REF+ALT (BAM)' | variable == 'Tumor REF (BAM)' | variable == 'Tumor ALT (BAM)'))), aes(x=factor(race),y=log(value+.001,base=2),fill=race)) + xlab("") + ylab("log2(Number of reads per patient)") + facet_grid(Cancer~variable) +
  geom_boxplot() +  labs(title=paste0("Exome read depth of *SNP* mutations in ", GENE,"\nfor patients WITH these ",GENE," mutations in MAF\n(mutations identified from MAF; read depth computed from BAMs)")) + 
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
  stat_summary(fun.data=n_fun,geom='text',size=3.5) + #calculate n
  ylim(-2.5,13) +
  stat_compare_means(size=5,method='t.test',mapping=aes(label=..p.signif..),
                     comparisons = list(c("Black","White")) ,
                     #label.y=1000, #adjust position of p-value
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("p<0.0001", "p<0.001", "p<0.01", "p<0.05", "NS")))

dev.off()

rm(list = ls())










#then the normal samples
pdf(paste0("/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/plots/",GENE,"_depth_from_MAF_by_race_normal.pdf"),height=14,width=13)

ggplot(subset(from_VCF_melted,(variable == 'Normal REF+ALT (MAF)' | variable == 'Normal REF (MAF)' | variable == 'Normal ALT (MAF)')), aes(x=factor(race),y=log(value+.001,base=2),fill=race)) + xlab("") + ylab("log2(Number of reads per patient)") + facet_grid(Cancer~variable) +
  geom_boxplot() +  labs(title=paste0("Exome read depth of mutations in ", GENE," for patients with these ",GENE," mutations in MAF")) + 
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
  #  ylim(-1.5,11) + 
  stat_summary(fun.data=n_fun,geom='text',size=3.5) + #calculate n
  ylim(-2.5,13) +
  stat_compare_means(size=5,method='t.test',mapping=aes(label=..p.signif..),
                     comparisons = list(c("Black","White")) ,
                     #label.y=1000, #adjust position of p-value
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("p<0.0001", "p<0.001", "p<0.01", "p<0.05", "NS")))

dev.off()
