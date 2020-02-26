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
rm(list = ls())

CANCER_LIST <- c('BRCA','LUAD','UCEC','KIRC','PRAD','COAD')
GENE <- 'TP53'

###############################################
#0. LOAD AND PROCESS DATA
############################################### 

#Retrieve gene coordinates from gtf file (downloaded from wget https://api.gdc.cancer.gov/data/25aa497c-e615-4cb7-8751-71f744f9691f; then paste <(zcat 25aa497c-e615-4cb7-8751-71f744f9691f | awk '$3 == "gene"' | cut -f1,4,5,7)  <(zcat 25aa497c-e615-4cb7-8751-71f744f9691f | awk '$3 == "gene"' | cut -d ' ' -f2,8)  > TCGA_GRCH38_v22.gtf
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

#Create dataframes to hold combined datasets
combined_data_master <- data.frame()
cdriver_master <- data.frame()
all_bams_all_mutations_master <- data.frame()

#Load data
for (CANCER_NAME in CANCER_LIST){
  #load overall depths data
  combined_data <- read.table(paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/overall_depths/',CANCER_NAME,'_overall_depths_etc.txt'),header=TRUE)
  
  #load gene depths data
  cdriver <- read.table(paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/gene_depths_pairs_from_VCF/',GENE,'/',CANCER_NAME,'_',names(GENE_AND_ENSEMBL_ID),'_depths_etc.txt'),header=TRUE)
  
  #load all bams all mutations in gene data
  all_bams_all_mutations <- read.table(paste0('/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/gene_depths_pairs_missing_from_VCF/',GENE,'/',CANCER_NAME,'_',names(GENE_AND_ENSEMBL_ID),'_depths_all_mutations_all_samples_etc.txt'),header=TRUE)
  all_bams_all_mutations$Cancer <- CANCER_NAME
  
  combined_data_master <- rbind(combined_data_master,combined_data)
  cdriver_master <- rbind(cdriver_master,cdriver)
  all_bams_all_mutations_master <- rbind(all_bams_all_mutations_master,all_bams_all_mutations)
}

#order by race
combined_data_master$race <- factor(combined_data_master$race, levels = c("White","Black","Asian","Unknown"),ordered=TRUE)  
cdriver_master$race <- factor(cdriver_master$race, levels = c("White","Black","Asian","Unknown")) 
all_bams_all_mutations_master$race <- factor(all_bams_all_mutations_master$race, levels = c("White","Black","Asian","Unknown")) 

#identify which from all_bams_all_mutations have at least 1 alt read
all_bams_all_mutations_absent_from_VCF_master <- all_bams_all_mutations_master[!interaction(all_bams_all_mutations_master[c('submitter_id','chromosome','mutation_position')]) %in% interaction(cdriver_master[c('submitter_id','chromosome','mutation_position')]),]

all_bams_all_mutations_absent_from_VCF_master <- subset(all_bams_all_mutations_absent_from_VCF_master,(alt_base_count_bam != 0))

#remove excess dataframes
rm(gene_subset)
rm(combined_data)
rm(cdriver)
rm(all_bams_all_mutations)


###############################################
#1. TOTAL NUMBER OF READS
############################################### 
pdf("/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/plots/total_exome_reads_per_patient.pdf",height=6,width=10)

n_fun <- function(x){return(data.frame(y=-12000000,label=paste(length(x),"patients")))} #function for calculating n

#EXOME
ggplot(combined_data_master, aes(x=factor(Cancer),y=read_depth_exome,fill=Cancer)) + xlab("") + ylab("Number of reads per patient") +
  geom_boxplot() +  labs(title=paste("Total number of exome reads per patient",sep='')) + 
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
  stat_summary(fun.data=n_fun,geom='text',size=5) # #calculate n 

dev.off()

###############################################
#2. TOTAL NUMBER OF EXOME READS BY RACE
############################################### 

#function to calculate n and place below each box
n_fun <- function(x){return(data.frame(y=-50000000,label=paste(length(x),"patients")))} 

#plot
pdf("/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/plots/total_exome_reads_per_patient_by_race.pdf",height=8,width=13)

ggplot(combined_data_master, aes(x=factor(race),y=read_depth_exome,fill=race)) + xlab("") + ylab("Number of reads per patient") + facet_wrap(~Cancer,ncol=4) + labs(title="Total number of exome reads per patient by race") + 
         geom_boxplot() + #fill=rep(c('#F8766D','#7CAE00','#00BFC4','#C77CFF'),length(CANCER_LIST))) +   
         theme(strip.text = element_text(size = 16,face='bold'),
              axis.text.x = element_blank(),
              axis.text.y = element_text(size=12),
              axis.title.y = element_text(size=14,face='bold',margin=margin(t=0,r=10,b=0,l=0)),
              axis.ticks= element_blank(), 
              plot.title = element_text(hjust=0.5,size=22,face='bold'),
              legend.position = 'bottom',legend.direction='horizontal',
              legend.text=element_text(size=18),
              legend.title=element_blank()) + 
              geom_point() + 
              ylim(-50000000,550000000) + 
         stat_summary(fun.data=n_fun,geom='text',size=3) + #calculate n
         scale_fill_manual(values=c('#F8766D','#7CAE00','#00BFC4','#C77CFF')) + #MANUALLY SET FILL COLORS
         
         #means comparisons
         stat_compare_means(size=4,method='t.test',mapping=aes(label=..p.signif..),
                            comparisons = list(c("Black","White"),c("Black","Asian"),c("White","Asian")),
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

pdf(paste0("/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/plots/",GENE,"_depth_from_VCF_by_race.pdf"),height=14,width=13)

#function for calculating n
n_fun <- function(x){return(data.frame(y=-2,label=paste(length(x),"mutation-\npatient pairs")))} 

#first, melt the dataframE
from_VCF_melted <- cdriver_master[,c("submitter_id","chromosome","mutation_position","race","Cancer","variant_read_depth_bam_MPILEUP","ref_base_count_bam","alt_base_count_bam")]
from_VCF_melted <- melt(from_VCF_melted,id.vars=c("submitter_id","chromosome","mutation_position","race","Cancer"),na.rm = TRUE)

#change variable names
from_VCF_melted$variable <- gsub("variant_read_depth_bam_MPILEUP","REF+ALT alleles",from_VCF_melted$variable)
from_VCF_melted$variable <- gsub("ref_base_count_bam","REF alleles",from_VCF_melted$variable)
from_VCF_melted$variable <- gsub("alt_base_count_bam","ALT alleles",from_VCF_melted$variable)

#reorder
from_VCF_melted$variable <- factor(from_VCF_melted$variable, levels = c("REF+ALT alleles","REF alleles","ALT alleles"),ordered=TRUE) 

ggplot(from_VCF_melted, aes(x=factor(race),y=log(value+1,base=2),fill=race)) + xlab("") + ylab("log2(Number of reads per patient)") + facet_grid(Cancer~variable) +
    geom_boxplot() +  labs(title=paste0("Exome read depth of mutations in ", GENE,"\nfor patients WITH these ",GENE," mutations in VCF\n(mutations identified from VCFs; read depth computed from BAMs)")) + 
  theme(strip.text = element_text(size = 14,face='bold'),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=14,face='bold',margin=margin(t=0,r=10,b=0,l=0)),
        axis.ticks= element_blank(), 
        plot.title = element_text(hjust=0.5,size=18,face='bold'),
        legend.position = 'bottom',legend.direction='horizontal',
        legend.text=element_text(size=14),
        legend.title=element_blank()) + 
        geom_point() + 
      #  ylim(-1.5,11) + 
    stat_summary(fun.data=n_fun,geom='text',size=3) + #calculate n
        ylim(-2.5,18) +
    stat_compare_means(size=3.3,method='t.test',mapping=aes(label=..p.signif..),
                       comparisons = list(c("Black","White"),c("Black","Asian"),c("White","Asian")) ,
                       #label.y=1000, #adjust position of p-value
                       symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("p<0.0001", "p<0.001", "p<0.01", "p<0.05", "NS")))
  
dev.off()


###############################################
#4. GENE READ-DEPTH BY RACE (calculated from bam) FOR PATIENTS *WITHOUT* MUTATIONS IN THAT GENE IN VCF
############################################### 

pdf(paste0("/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/plots/",GENE,"_depth_absent_from_VCF_by_race.pdf"),height=14,width=13)

#function for calculating n
n_fun <- function(x){return(data.frame(y=-2,label=paste(length(x),"mutation-\npatient pairs")))} 

#first, melt the dataframe....
missing_from_VCF_melted <- all_bams_all_mutations_absent_from_VCF_master[,c("submitter_id","chromosome","mutation_position","race","Cancer","variant_read_depth_bam_MPILEUP","ref_base_count_bam","alt_base_count_bam")]

missing_from_VCF_melted <- melt(missing_from_VCF_melted,id.vars=c("submitter_id","chromosome","mutation_position","race","Cancer"),na.rm = TRUE)

#change variable names
missing_from_VCF_melted$variable <- gsub("variant_read_depth_bam_MPILEUP","REF+ALT alleles",missing_from_VCF_melted$variable)
missing_from_VCF_melted$variable <- gsub("ref_base_count_bam","REF alleles",missing_from_VCF_melted$variable)
missing_from_VCF_melted$variable <- gsub("alt_base_count_bam","ALT alleles",missing_from_VCF_melted$variable)

#reorder
missing_from_VCF_melted$variable <- factor(missing_from_VCF_melted$variable, levels = c("REF+ALT alleles","REF alleles","ALT alleles"),ordered=TRUE) 
missing_from_VCF_melted$Cancer <- factor(missing_from_VCF_melted$Cancer, levels = CANCER_LIST,ordered=TRUE) 

ggplot(missing_from_VCF_melted, aes(x=factor(race),y=log(value+1,base=2),fill=race)) + xlab("") + ylab("log2(Number of reads per patient)") + facet_grid(Cancer~variable) +
  geom_boxplot() +  labs(title=paste0("Exome read depth of mutations in ", GENE,"\nfor patients WITHOUT these ",GENE," mutations in VCF\n(mutations identified from VCFs; read depth computed from BAMs)")) + 
  theme(strip.text = element_text(size = 14,face='bold'),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=14,face='bold',margin=margin(t=0,r=10,b=0,l=0)),
        axis.ticks= element_blank(), 
        plot.title = element_text(hjust=0.5,size=18,face='bold'),
        legend.position = 'bottom',legend.direction='horizontal',
        legend.text=element_text(size=14),
        legend.title=element_blank()) + 
  geom_point() + 
  #  ylim(-1.5,11) + 
  stat_summary(fun.data=n_fun,geom='text',size=3) + #calculate n
  ylim(-2.5,18) +
  stat_compare_means(size=3.3,method='t.test',mapping=aes(label=..p.signif..),
                     comparisons = list(c("Black","White"),c("Black","Asian"),c("White","Asian")) ,
                     #label.y=1000, #adjust position of p-value
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("p<0.0001", "p<0.001", "p<0.01", "p<0.05", "NS")))

dev.off()

  
###############################################
#5. ALT ALLELE CONCENTRATION FOR PATIENTS *WITH* MUTATIONS IN THAT GENE IN VCF
############################################### 
pdf(paste0("/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/plots/",GENE,"_alt_allele_concentration_pairs_from_VCF.pdf"),height=9,width=13)
  
#function for calculating n
n_fun <- function(x){return(data.frame(y=-0.09,label=paste(length(x),"mutation-\npatient pairs")))} 

ggplot(cdriver_master, aes(x=factor(race),y=alt_allele_concent_bam,fill=race)) + xlab("") + ylab("Alt. allele concentration") + facet_wrap(~Cancer) +
  geom_boxplot() +  labs(title=paste0("Alt allele concentration in ", GENE,"\nfor patients WITH these ", GENE," mutations in VCF\n(mutations identified from VCF; alt allele concentration computed from BAM)")) + 
  theme(strip.text = element_text(size = 16,face='bold'),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=14,face='bold',margin=margin(t=0,r=10,b=0,l=0)),
        axis.ticks= element_blank(), 
        plot.title = element_text(hjust=0.5,size=18,face='bold'),
        legend.position = 'bottom',legend.direction='horizontal',
        legend.text=element_text(size=18),
        legend.title=element_blank()) + 
  geom_point() + 
  #  ylim(-1.5,11) + 
  stat_summary(fun.data=n_fun,geom='text',size=3) + #calculate n
  ylim(-0.10,1.3) +
  stat_compare_means(size=4,method='t.test',mapping=aes(label=..p.signif..),
                     comparisons = list(c("Black","White"),c("Black","Asian"),c("White","Asian")) ,
                     #label.y=1000, #adjust position of p-value
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("p<0.0001", "p<0.001", "p<0.01", "p<0.05", "NS")))

dev.off()  
  
  
###############################################
#5. ALT ALLELE CONCENTRATION FOR PATIENTS *WITHOUT* THOSE MUTATIONS IN THAT GENE IN VCF
############################################### 

#why are some larger than 1????

pdf(paste0("/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/plots/",GENE,"_alt_allele_concentration_pairs_absent_from_VCF.pdf"),height=9,width=13)

#function for calculating n
n_fun <- function(x){return(data.frame(y=-0.09,label=paste(length(x),"mutation-\npatient pairs")))} 

ordered <- all_bams_all_mutations_absent_from_VCF_master
ordered$Cancer <- factor(ordered$Cancer, levels = CANCER_LIST,ordered=TRUE)  

ggplot(subset(all_bams_all_mutations_absent_from_VCF_master,alt_allele_concent_bam <= 1.0), aes(x=factor(race),y=alt_allele_concent_bam,fill=race)) + xlab("") + ylab("Alt. allele concentration") + facet_wrap(~Cancer) +
  geom_boxplot() +  labs(title=paste0("Alt allele concentration in ", GENE,"\nfor patients WITHOUT these ", GENE," mutations in VCF\n(mutations identified from VCF; alt allele concentration computed from BAM)")) + 
  theme(strip.text = element_text(size = 16,face='bold'),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=14,face='bold',margin=margin(t=0,r=10,b=0,l=0)),
        axis.ticks= element_blank(), 
        plot.title = element_text(hjust=0.5,size=18,face='bold'),
        legend.position = 'bottom',legend.direction='horizontal',
        legend.text=element_text(size=18),
        legend.title=element_blank()) + 
  geom_point() + 
  #  ylim(-1.5,11) + 
  stat_summary(fun.data=n_fun,geom='text',size=3) + #calculate n
  ylim(-0.10,1.3) +
  stat_compare_means(size=4,method='t.test',mapping=aes(label=..p.signif..),
                     comparisons = list(c("Black","White"),c("Black","Asian"),c("White","Asian")) ,
                     #label.y=1000, #adjust position of p-value
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("p<0.0001", "p<0.001", "p<0.01", "p<0.05", "NS")))

dev.off()  
  
  
  
###############################################
#6. TOTAL NUMBER OF READS PER PATIENT BY SOURCE SITE AND RACE
############################################### 
  
  #function for calculating n
  n_fun <- function(x){return(data.frame(y=-9000,label=paste(length(x),"patients")))} 
  
for (CANCER_NAME in CANCER_LIST){  
  
  pdf(paste0("/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/plots/overall_depth_by_source_site_",CANCER_NAME,".pdf"),height=9,width=11)
  
    source_site <- read.table(paste0("/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/overall_depths/",CANCER_NAME,"_overall_depths_etc.txt"),header=TRUE)
    
    #subset only the white and black samples
    source_site <- subset(source_site, (race=='White') | (race=='Black'))
    
    #keep only the centers with at least 10 samples
    source_site <- source_site[source_site$Source_Site_exome %in% names(which(table(source_site$Source_Site_exome) >= 10)),]
    
    #change variable names
    source_site$Source_Site_exome <- gsub("_"," ",source_site$Source_Site_exome)
    
   TSS <- ggplot(source_site, aes(fill=race,x=race, y=read_depth_exome)) + # use 'fill' to set colors
      geom_boxplot(outlier.shape=NA)+  ylab("Exome read depth") +  xlab("") + labs(title = paste0("CANCER: ",CANCER_NAME,"\nTotal number of exome reads per patient from source sites contributing at least 10 samples")) +  
      theme(plot.title = element_text(hjust=0.5),
            axis.title.y = element_text(size=14,margin = margin(t = 0, r = 5, b = 0, l = 0)),
            axis.title.x = element_text(margin = margin(t = 15, r = 10, b = 5, l = 10)),
            axis.text.x = element_blank(),axis.ticks= element_blank(),
            axis.text.y = element_text(size=12),
            legend.text=element_text(size=12), legend.title=element_blank(),
            legend.background = element_rect(linetype='solid', colour='black'),
            legend.position = c(0.5,-0.04),legend.direction='horizontal',plot.margin = margin(t=5,r=15,b=15,l=5))+
      geom_point() + stat_summary(fun.data=n_fun,geom='text',size=4) + #calculate n
      ylim(-10000,max(source_site$read_depth_exome)+50000000) +
      facet_wrap(~Source_Site_exome, strip.position="top",ncol=5) + theme(strip.text = element_text(face="bold", size=10)) +
      stat_compare_means(size=3.3,method='t.test',mapping=aes(label=..p.signif..),
                         comparisons = list(c("Black","White")),
                         
                         symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("p<0.0001", "p<0.001", "p<0.01", "p<0.05", "NS"))) 
    print(TSS)  
    dev.off()
    
  }  
  
###############################################
#7. TOTAL NUMBER OF READS PER PATIENT BY SEQUENCING CENTER
############################################### 
  
pdf(paste0("/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/plots/overall_depth_by_seq_center.pdf"),height=13,width=13)

#function to calculate n and place below each box
n_fun <- function(x){return(data.frame(y=0,label=paste(length(x),"patients")))} 

#change variable names (center codes from https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/center-codes)
combined_data_master$center_exome <- gsub("10","Baylor",combined_data_master$center_exome)
combined_data_master$center_exome <- gsub("9","WashU",combined_data_master$center_exome)
combined_data_master$center_exome <- gsub("8","Broad",combined_data_master$center_exome)
combined_data_master$center_exome <- gsub("1","Broad",combined_data_master$center_exome)


ggplot(subset(combined_data_master, (center_exome == 'Broad' | center_exome == 'WashU' | center_exome == 'Baylor')), aes(x=factor(race),y=read_depth_exome,fill=race)) + xlab("") + ylab("Number of reads") + facet_grid(Cancer~center_exome) +
  geom_boxplot() +  labs(title="Total number of exome reads per patient by sequencing center") + 
  theme(strip.text = element_text(size = 16,face='bold'),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=14,face='bold',margin=margin(t=0,r=10,b=0,l=0)),
        axis.ticks= element_blank(), 
        plot.title = element_text(hjust=0.5,size=18,face='bold'),
        legend.position = 'bottom',legend.direction='horizontal',
        legend.text=element_text(size=18),
        legend.title=element_blank()) + 
  geom_point() + 
  ylim(-1.5,550000000) + 
  stat_summary(fun.data=n_fun,geom='text',size=3) + #calculate n
  #ylim(-0.10,1.3) +
  stat_compare_means(size=4,method='t.test',mapping=aes(label=..p.signif..),
                     comparisons = list(c("Black","White"),c("Black","Asian"),c("White","Asian")) ,
                     #label.y=1000, #adjust position of p-value
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("p<0.0001", "p<0.001", "p<0.01", "p<0.05", "NS")))

dev.off()    








###############################################
#8. GENE READ-QUALITY BY RACE (calculated from bam) FOR PATIENTS *WITH* MUTATIONS IN THAT GENE IN VCF
############################################### 

pdf(paste0("/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/plots/",GENE,"_quality_from_VCF_by_race.pdf"),height=14,width=9)

#function for calculating n
n_fun <- function(x){return(data.frame(y=-6,label=paste(length(x),"mutation-\npatient pairs")))} 

#first, melt the dataframE
from_VCF_melted <- cdriver_master[,c("submitter_id","chromosome","mutation_position","race","Cancer","avg_ref_base_quality_bam","avg_alt_base_quality_bam")]
from_VCF_melted <- melt(from_VCF_melted,id.vars=c("submitter_id","chromosome","mutation_position","race","Cancer"),na.rm = TRUE)

#change variable names
from_VCF_melted$variable <- gsub("avg_ref_base_quality_bam","REF alleles",from_VCF_melted$variable)
from_VCF_melted$variable <- gsub("avg_alt_base_quality_bam","ALT alleles",from_VCF_melted$variable)

#reorder
from_VCF_melted$variable <- factor(from_VCF_melted$variable, levels = c("REF alleles","ALT alleles"),ordered=TRUE) 

ggplot(from_VCF_melted, aes(x=factor(race),y=value,fill=race)) + xlab("") + ylab("Base quality per patient") + facet_grid(Cancer~variable) +
  geom_boxplot() +  labs(title=paste0("Average base quality of mutations in ", GENE,"\nfor patients WITH these ",GENE," mutations in VCF\n(mutations identified from VCFs; base quality from BAMs)")) + 
  theme(strip.text = element_text(size = 14,face='bold'),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=14,face='bold',margin=margin(t=0,r=10,b=0,l=0)),
        axis.ticks= element_blank(), 
        plot.title = element_text(hjust=0.5,size=18,face='bold'),
        legend.position = 'bottom',legend.direction='horizontal',
        legend.text=element_text(size=14),
        legend.title=element_blank()) + 
  geom_point() + 
  stat_summary(fun.data=n_fun,geom='text',size=3) + #calculate n
  ylim(-8,60) +
  stat_compare_means(size=3.3,method='t.test',mapping=aes(label=..p.signif..),
                     comparisons = list(c("Black","White"),c("Black","Asian"),c("White","Asian")) ,
                     #label.y=1000, #adjust position of p-value
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("p<0.0001", "p<0.001", "p<0.01", "p<0.05", "NS")))



dev.off()



###############################################
#9. GENE READ-QUALITY BY RACE (calculated from bam) FOR PATIENTS *WITHOUT* MUTATIONS IN THAT GENE IN VCF
############################################### 

pdf(paste0("/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/plots/",GENE,"_quality_absent_from_VCF_by_race.pdf"),height=14,width=9)

#function for calculating n
n_fun <- function(x){return(data.frame(y=-14,label=paste(length(x),"mutation-\npatient pairs")))} 

#first, melt the dataframE
from_VCF_melted <- all_bams_all_mutations_absent_from_VCF_master[,c("submitter_id","chromosome","mutation_position","race","Cancer","avg_ref_base_quality_bam","avg_alt_base_quality_bam")]
from_VCF_melted <- melt(from_VCF_melted,id.vars=c("submitter_id","chromosome","mutation_position","race","Cancer"),na.rm = TRUE)

#change variable names
from_VCF_melted$variable <- gsub("avg_ref_base_quality_bam","REF alleles",from_VCF_melted$variable)
from_VCF_melted$variable <- gsub("avg_alt_base_quality_bam","ALT alleles",from_VCF_melted$variable)

#reorder
from_VCF_melted$variable <- factor(from_VCF_melted$variable, levels = c("REF alleles","ALT alleles"),ordered=TRUE) 

ggplot(from_VCF_melted, aes(x=factor(race),y=value,fill=race)) + xlab("") + ylab("Base quality per patient") + facet_grid(Cancer~variable) +
  geom_boxplot() +  labs(title=paste0("Average base quality of mutations in ", GENE,"\nfor patients WITHOUT these ",GENE," mutations in VCF\n(mutations identified from VCFs; base quality from BAMs)")) + 
  theme(strip.text = element_text(size = 14,face='bold'),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=14,face='bold',margin=margin(t=0,r=10,b=0,l=0)),
        axis.ticks= element_blank(), 
        plot.title = element_text(hjust=0.5,size=18,face='bold'),
        legend.position = 'bottom',legend.direction='horizontal',
        legend.text=element_text(size=14),
        legend.title=element_blank()) + 
  geom_point() + 
  stat_summary(fun.data=n_fun,geom='text',size=3) + #calculate n
  ylim(-17,105) +
  stat_compare_means(size=3.3,method='t.test',mapping=aes(label=..p.signif..),
                     comparisons = list(c("Black","White"),c("Black","Asian"),c("White","Asian")) ,
                     #label.y=1000, #adjust position of p-value
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("p<0.0001", "p<0.001", "p<0.01", "p<0.05", "NS")))


dev.off()





  
  
  
  

  
  
  
  
  
  
  
  




###############################################
#NUMBER OF SOMATIC MUTATIONS PER PATIENT
############################################### 
n_fun <- function(x){return(data.frame(y=-0.5,label=paste(length(x),"patients")))} #function for calculating n

system(paste0("mkdir /home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/plots/mutations_per_patient/",GENE))

for (CANCER_NAME in CANCER_LIST){  

  x <- subset(cdriver_master,(Cancer==CANCER_NAME))
  
  #comparisons depend on which races present
  if ("Asian" %in% as.character(unique(x$race)) & "Black" %in% as.character(unique(x$race)) & "White" %in% as.character(unique(x$race))){
    comparisons_cmnd <- list(c("Black","White"),c("Black","Asian"),c("White","Asian"))  
  }
  if (!("Asian" %in% as.character(unique(x$race))) & "Black" %in% as.character(unique(x$race)) & "White" %in% as.character(unique(x$race))){
    comparisons_cmnd <- list(c("Black","White")) 
  }
  if ("Asian" %in% as.character(unique(x$race)) & !("Black" %in% as.character(unique(x$race))) & "White" %in% as.character(unique(x$race))){
    comparisons_cmnd <- list(c("Black","White"),c("Black","Asian"),c("White","Asian"))
  }
  
  assign(paste0('mutations_per_patient_',CANCER_NAME), ggplot(x[!duplicated(x$submitter_id),], aes(x=factor(race),y=as.numeric(mutations_count_for_caseID),fill=race)) + xlab("") + ylab("Number of mutatant genotypes") +
    geom_violin() +  labs(title=CANCER_NAME) + 
    theme(strip.text = element_text(size = 8),axis.text.x = element_blank(),axis.ticks= element_blank(), plot.title = element_text(hjust=0.5,size=16), legend.text=element_text(size=12),legend.title=element_blank()) + geom_point() + 
    stat_summary(fun.data=n_fun,geom='text',size=4) + #calculate n
    guides(fill = guide_legend(reverse=FALSE)) + #reverse order of legend
    
    stat_compare_means(size=3.3,method='t.test',mapping=aes(label=..p.signif..),
                       comparisons = comparisons_cmnd,
                       #label.y=5, #adjust position of p-value
                       symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("p<0.0001", "p<0.001", "p<0.01", "p<0.05", "NS"))))
}

all_plots <- ggarrange(mutations_per_patient_BRCA,mutations_per_patient_LUAD,mutations_per_patient_LUSC,mutations_per_patient_THCA,common.legend=TRUE,legend='bottom',nrow=2,ncol=2)

pdf(paste0("/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/plots/mutations_per_patient/",GENE,"/",GENE,"_mutations_per_patient.pdf"),height=9,width=10)
print(annotate_figure(all_plots,top=text_grob(paste("Number of ",sub('^([^_]*).*', '\\1',names(DRIVER_GENES)[1])," mutatant genotypes per patient \n(in patients with at least 1 ",sub('^([^_]*).*', '\\1',names(DRIVER_GENES)[1])," mutant genotype)",sep=''), color = "blue", face = "bold", size = 18)))
dev.off()














#############################################################

#9. MUTATION POSITION

#LINEAR PLOT OF COUNTS ALONG GENE BY RACE
positions <- table(cdriver_master$race,cdriver_master$mutation_position)

#get list of known variants
known_pathogenic_variants <- scan('/Users/m187735/Desktop/Rscripts/p53_pathogenic.txt', character(), quote = "")
detected_variants <- unique(cdriver_master$mutation_position)
overlap <- intersect(known_pathogenic_variants,detected_variants)    
overlapping_position_index <- match(overlap,colnames(positions))
non_overlap <- detected_variants [ ! detected_variants %in% known_pathogenic_variants ]
nonoverlapping_position_index <- match(non_overlap,colnames(positions))

linear_plot <- barplot(positions,las=2,col=c("red","green","blue"),cex.names=0.25, ylab= "Number of mutations",xlab="", legend = rownames(positions),yaxt='n',xaxt='n',main=paste("Cancer: ",CANCER_NAME,"\nCounts for ",sub('^([^_]*).*', '\\1',names(DRIVER_GENES)[1])," mutations\n(source: somatic VCFs)",sep=''),cex.main=1)
axis(2, at=seq(0,34,2))
axis(1, at=linear_plot[overlapping_position_index],labels = colnames(positions)[overlapping_position_index],col.axis='orange',las=2,cex.axis=0.85)
axis(1, at=linear_plot[nonoverlapping_position_index],labels = colnames(positions)[nonoverlapping_position_index],col.axis='black',las=2,cex.axis=0.85)


#8. LINEAR PLOT OF TOTAL READ DEPTH

#from BAM (in patients with TP53 mutations called in VCF)
positions <- cdriver_master[,c(3,16,35)]
positions <- aggregate(alt_base_count_bam~mutation_position+race,data=positions,FUN=sum)
positions <- xtabs(alt_base_count_bam ~ mutation_position + race,data=positions)

linear_plots <- barplot(t(positions),las=2,col=c("red","green","blue"),cex.names=0.25, ylab= "Total number of ALT bases at position",xlab="", legend = rownames(t(positions)),yaxt='n',xaxt='n',main=paste("Cancer: ",CANCER_NAME,"\nTotal read depth of ALT alleles of ",length(unique(cdriver_master$mutation_position))," mutated loci in ",sub('^([^_]*).*', '\\1',names(DRIVER_GENES)[1])," \nfor locus-patient pairs with mutations called in VCF\n(source: tumor exome BAMs)",sep=''),cex.main=1)
axis(2, at=seq(0,120000,100))
axis(1, at=linear_plots[overlapping_position_index],labels = colnames(t(positions))[overlapping_position_index],col.axis='orange',las=2,cex.axis=0.85)
axis(1, at=linear_plots[nonoverlapping_position_index],labels = colnames(t(positions))[nonoverlapping_position_index],col.axis='black',las=2,cex.axis=0.85)    


#from BAM (in ALL patients, regardless of their VCF mutations)
positions <- cdriver_master_all_bams_all_mutations[,c(3,12,20)]
positions <- aggregate(alt_base_count_bam~mutation_position+race,data=positions,FUN=sum)
positions <- xtabs(alt_base_count_bam ~ mutation_position + race,data=positions)

linear_plots <- barplot(t(positions),las=2,col=c("red","green","blue"),cex.names=0.25, ylab= "Total number of ALT bases at position",xlab="", legend = rownames(t(positions)),yaxt='n',xaxt='n',main=paste("Cancer: ",CANCER_NAME,"\nTotal read depth of ALT alleles of ",length(unique(cdriver_master$mutation_position))," mutated loci in ",sub('^([^_]*).*', '\\1',names(DRIVER_GENES)[1])," \nfor ALL locus-patient pairs (regardless of whether mutations called in VCF)\n(source: tumor exome BAMs)",sep=''),cex.main=1)
axis(2, at=seq(0,120000,100))
axis(1, at=linear_plots[overlapping_position_index],labels = colnames(t(positions))[overlapping_position_index],col.axis='orange',las=2,cex.axis=0.85)
axis(1, at=linear_plots[nonoverlapping_position_index],labels = colnames(t(positions))[nonoverlapping_position_index],col.axis='black',las=2,cex.axis=0.85)    

#from BAM (for locus-patient pairs missing from VCF)
positions <- driver_all_bams_all_mutations_minus_VCF_mutated[,c(2,6,12)]
positions <- aggregate(alt_base_count_bam~mutation_position+race,data=positions,FUN=sum)
positions <- xtabs(alt_base_count_bam ~ mutation_position + race,data=positions)

linear_plots <- barplot(t(positions),las=2,col=c("red","green","blue"),cex.names=0.25, ylab= "Total number of ALT bases at position",xlab="", legend = rownames(t(positions)),yaxt='n',xaxt='n',main=paste("Cancer: ",CANCER_NAME,"\nTotal read depth of ALT alleles of ",length(unique(cdriver_master$mutation_position))," mutated loci in ",sub('^([^_]*).*', '\\1',names(DRIVER_GENES)[1])," \nfor locus-patient pairs without mutation called in VCF\n(source: tumor exome BAMs)",sep=''),cex.main=1)
axis(2, at=seq(0,120000,100))
axis(1, at=linear_plots[overlapping_position_index],labels = colnames(t(positions))[overlapping_position_index],col.axis='orange',las=2,cex.axis=0.85)
axis(1, at=linear_plots[nonoverlapping_position_index],labels = colnames(t(positions))[nonoverlapping_position_index],col.axis='black',las=2,cex.axis=0.85)    



