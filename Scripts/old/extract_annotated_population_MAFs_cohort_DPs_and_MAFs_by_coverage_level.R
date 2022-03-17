library(stringr)
library(reshape2)
library(iterators)
library(plyr)
library(dplyr)

args<-commandArgs(TRUE)
.libPaths()
if (length(args)!=3) {
  stop("Usage example: Rscript --vanilla all_proteins_all_chr/BRCA/chr3 Black Primary_Tumor", call.=FALSE)
} else if (length(args)==3) {
  # default output file
  PROJECT <- toString(args[1])
  RACE <- toString(args[2])
  TISSUE <- toString(args[3])
}



filename=paste0("/research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/processing/TCGA/processed_data/annotate_PASS_TCGA/",PROJECT,"/",TISSUE,"_",RACE,"_final_VQSR_chr_PASS.vcf.hg38_multianno.vcf.recode.vcf") #remove reoode if not doing both races together
N_ROWS <- system(paste("wc -l ",filename," | awk '{print $1}'"),intern=TRUE)

OUTPUT_FILE <- paste0("/research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/processing/TCGA/processed_data/annotate_PASS_TCGA/",PROJECT,"/",TISSUE,"_",RACE,"_final_VQSR_PASS_population_MAFs_cohort_DP_and_MAF_by_coverage_level.txt")
#system(paste0("rm ", OUTPUT_FILE))

#read file line by line
connection <- ireadLines(filename,n=1)
ALL_ROWS <- data.frame()

for (i in (1:N_ROWS)){
  print(i)
  ROW <- data.frame(nextElem(connection,))
  
  
  
  if (grepl('#', ROW$nextElem.connection..., fixed = TRUE) == TRUE){
    sample_IDs <- ROW
    
    try(sample_IDs <- ldply(strsplit(as.character(sample_IDs$nextElem.connection...), split = "\t")),silent = TRUE)
    try(sample_IDs <- sample_IDs[1,10:ncol(sample_IDs)],silent=TRUE)
    
    rm(ROW)
    i <- i+1
    next
  }
  
  
  ROW <- ldply(strsplit(as.character(ROW$nextElem.connection.), split = "\t"))
  names(ROW)[10:ncol(ROW)] <- sample_IDs
  ROW_GT <- ROW
  
  #isolate AD, GT from each column
  ROW[,c(10:ncol(ROW))] <- sapply(ROW[,c(10:ncol(ROW))], function (x) sub("^([^:]*:[^:]*:[^:]*).*", "\\1",x))
  ROW[,c(10:ncol(ROW))] <- sapply(ROW[,c(10:ncol(ROW))], function (x) gsub("[:|,]","|",x)) #replace punctuation with pipe
  # ROW[,c(10:ncol(ROW))] <- sapply(ROW[,c(10:ncol(ROW))], function (x) sub(".*,.*", "0",x)) #if any samples had no reads
  #ROW[,c(10:ncol(ROW))] <- sapply(ROW[,c(10:ncol(ROW))], function (x) sub("\\./\\.", "0",x)) #if any samples had no reads
  # ROW[,c(10:ncol(ROW))] <- sapply(ROW[,c(10:ncol(ROW))], function (x) as.numeric(x)) #must be numeric
  
  
  
  #change FORMAT and rm extraneous cols
  ROW[,9] <- 'GT|REF|ALT|DP'
  ROW <- ROW[,c(1,2,4,5,8,9,10:(ncol(ROW)))]
  names(ROW)[1:6] <- c('CHR','POS','REF','ALT','INFO','FORMAT')
  
  #isolate VQSLOD
  ROW$VQSLOD <- str_extract(ROW$INFO,"(?<=VQSLOD=)[^;]*(?=;|$)")
  
  #isolate gene name
  ROW$Gene <- str_extract(ROW$INFO, "(?<=Gene.refGene=)[^;]*(?=;|$)")
  
  #isolate MQ
  ROW$MQ <- str_extract(ROW$INFO, "(?<=MQ=)[^;]*(?=;|$)")
  
  #isolate population allele freqs of interest and rearrange
  ROW$EUR_1kGenomes <- str_extract(ROW$INFO, "(?<=EUR.sites.2015_08=)[^;]*(?=;|$)")
  ROW$AFR_1kGenomes <- str_extract(ROW$INFO, "(?<=AFR.sites.2015_08=)[^;]*(?=;|$)")
  
  ROW <- ROW[,c(1:4,6,(ncol(ROW)-4):ncol(ROW),7:(ncol(ROW)-5))]
  
  #create dataframe of samples only
  ALL_PATIENTS <- ROW[,c(11:ncol(ROW))]
  
  #calculate number of samples (for race) that are homo ref, hetero, homo alt, missing
  ROW$count_all_patients <- length(ALL_PATIENTS)
  ROW$count_all_REF_reads <- sum(as.numeric(sapply(ALL_PATIENTS, function (x) (strsplit(x,split='\\|')[[1]][2]))))
  ROW$count_all_ALT_reads <- sum(as.numeric(sapply(ALL_PATIENTS, function (x) (strsplit(x,split='\\|')[[1]][3]))))
  ROW$count_all_reads <- sum(as.numeric(sapply(ALL_PATIENTS, function (x) (strsplit(x,split='\\|')[[1]][4]))))
  
  ROW$GT_.._from_all_patients <- sum(grepl("\\./\\.",sapply(ALL_PATIENTS, function (x) (strsplit(x,split='\\|')[[1]][1]))))
  ROW$GT_00_from_all_patients <- sum(grepl("0/0",sapply(ALL_PATIENTS, function (x) (strsplit(x,split='\\|')[[1]][1]))))
  ROW$GT_01_from_all_patients <- sum(grepl("0/1",sapply(ALL_PATIENTS, function (x) (strsplit(x,split='\\|')[[1]][1]))))
  ROW$GT_11_from_all_patients <- sum(grepl("1/1",sapply(ALL_PATIENTS, function (x) (strsplit(x,split='\\|')[[1]][1]))))
  
  #calculate BAF
  ROW$BAF_all_patients <- as.numeric(2*ROW$GT_11_from_all_patients + ROW$GT_01_from_all_patients) / as.numeric(2*ROW$GT_11_from_all_patients + 2*ROW$GT_01_from_all_patients + 2*ROW$GT_00_from_all_patients)
  ROW$BAF_all_patients_incl_.. <- as.numeric(2*ROW$GT_11_from_all_patients + ROW$GT_01_from_all_patients) / as.numeric(2*ROW$GT_11_from_all_patients + 2*ROW$GT_01_from_all_patients + 2*ROW$GT_00_from_all_patients + 2*ROW$GT_.._from_all_patients)
  
  
  #isolate patients, getting rid of any ./. genotypes
  # PATIENTS <- ALL_PATIENTS[sapply(ALL_PATIENTS, function (x) (strsplit(x,split='\\|')[[1]][1] != './.'))] 
  
  PATIENTS <- ALL_PATIENTS
  
  PATIENTS_placeholder <- sapply(PATIENTS, function (x) (as.numeric(strsplit(x,split='\\|')[[1]][4] ) == "placeholder")) #remove any ./. in DP col
      
  PATIENTS <- PATIENTS[,sapply(PATIENTS_placeholder, function (x) (!is.na(x)))]
  
  #  if (length(PATIENTS) == 1){
  #    PATIENTS <- PATIENTS[sapply(PATIENTS, function (x) (strsplit(x,split='\\|')[[1]][4] != '.'))] #needs to have some alt reads
  #  }
  #  if (length(PATIENTS) > 1){
  #    PATIENTS <- PATIENTS[,sapply(PATIENTS, function (x) (strsplit(x,split='\\|')[[1]][4] != '.'))] #needs to have some alt reads
  #  }
  
  #only interested in the loci that are biallelic
  if(length(PATIENTS) > 1 & grepl(",",ROW$REF) == FALSE & grepl(",",ROW$ALT) == FALSE){
    
    #  PATIENTS <- PATIENTS[,sapply(PATIENTS, function (x) (as.numeric(strsplit(x,split='\\|')[[1]][4] ) <= DP))]
    
    
    #calcualte the same thing for individuals with < xx reads at this position
    #identify patients with alternative-supporting reads at each position
    
    BAF_BY_COVERAGE <- function(DP1,DP2,GREATER_OR_LESS_THAN){
      
      if (DP1==DP2){
        if (GREATER_OR_LESS_THAN == '<='){
          COV_INDIVIDUALS <- PATIENTS[,sapply(PATIENTS, function (x) (as.numeric(strsplit(x,split='\\|')[[1]][4] ) <= DP1))] #samples with DP below
          DIRECTION <- 'below'
        }  
        
        if (GREATER_OR_LESS_THAN == '>='){
          COV_INDIVIDUALS <- PATIENTS[,sapply(PATIENTS, function (x) (as.numeric(strsplit(x,split='\\|')[[1]][4] ) >= DP1))] #samples with DP above
          DIRECTION <- 'above'
        }  
      }
      
      if (DP1 != DP2){
        if (GREATER_OR_LESS_THAN == '<='){
          COV_INDIVIDUALS <- PATIENTS[,sapply(PATIENTS, function (x) (as.numeric(strsplit(x,split='\\|')[[1]][4] ) >= DP1))] 
          if(length(COV_INDIVIDUALS)>1){COV_INDIVIDUALS <- COV_INDIVIDUALS[,sapply(COV_INDIVIDUALS, function (x) (as.numeric(strsplit(x,split='\\|')[[1]][4] ) <= DP2))]}
          if(length(COV_INDIVIDUALS)==1){COV_INDIVIDUALS <- COV_INDIVIDUALS[as.numeric(strsplit(COV_INDIVIDUALS,split='\\|')[[1]][4] ) <= DP2]}
          DIRECTION <- 'within'
        }  
      }
      
      
      
      count_patients_with_DP_ <- length(COV_INDIVIDUALS)
      
      #count genotypes
      GT_00_patients_with_DP_ <- sum(grepl("0/0",sapply(COV_INDIVIDUALS, function (x) (strsplit(x,split='\\|')[[1]][1]))))
      GT_01_patients_with_DP_ <- sum(grepl("0/1",sapply(COV_INDIVIDUALS, function (x) (strsplit(x,split='\\|')[[1]][1]))))
      GT_11_patients_with_DP_ <- sum(grepl("1/1",sapply(COV_INDIVIDUALS, function (x) (strsplit(x,split='\\|')[[1]][1]))))
      GT_.._patients_with_DP_ <- sum(grepl("\\./\\.",sapply(COV_INDIVIDUALS, function (x) (strsplit(x,split='\\|')[[1]][1]))))
      
      #calculate genotype proportions
      #  PROP_GT_00_patients_with_DP_ <- GT_00_patients_with_DP_ / (GT_00_patients_with_DP_ + GT_01_patients_with_DP_ + GT_11_patients_with_DP_ + GT_.._patients_with_DP_)
      #  PROP_GT_01_patients_with_DP_ <- GT_01_patients_with_DP_ / (GT_00_patients_with_DP_ + GT_01_patients_with_DP_ + GT_11_patients_with_DP_ + GT_.._patients_with_DP_)
      #  PROP_GT_11_patients_with_DP_ <- GT_11_patients_with_DP_ / (GT_00_patients_with_DP_ + GT_01_patients_with_DP_ + GT_11_patients_with_DP_ + GT_.._patients_with_DP_)
      #  PROP_GT_.._patients_with_DP_ <- GT_.._patients_with_DP_ / (GT_00_patients_with_DP_ + GT_01_patients_with_DP_ + GT_11_patients_with_DP_ + GT_.._patients_with_DP_)
      
      #calculate BAF
      BAF_patients_with_DP_ <- as.numeric(2*GT_11_patients_with_DP_ + GT_01_patients_with_DP_) / as.numeric(2*GT_11_patients_with_DP_ + 2*GT_01_patients_with_DP_ + 2*GT_00_patients_with_DP_)
      BAF_patients_incl_.._with_DP_ <- as.numeric(2*GT_11_patients_with_DP_ + GT_01_patients_with_DP_) / as.numeric(2*GT_11_patients_with_DP_ + 2*GT_01_patients_with_DP_ + 2*GT_00_patients_with_DP_ + 2*GT_.._patients_with_DP_)
      
      #add these data to ROW
      #  ROW <- cbind(ROW, count_patients_with_DP_,GT_00_patients_with_DP_,PROP_GT_00_patients_with_DP_,GT_01_patients_with_DP_,PROP_GT_01_patients_with_DP_,GT_11_patients_with_DP_,PROP_GT_11_patients_with_DP_,GT_.._patients_with_DP_,PROP_GT_.._patients_with_DP_,BAF_patients_with_DP_,BAF_patients_incl_.._with_DP_)
      
      ROW <- cbind(ROW, count_patients_with_DP_,GT_00_patients_with_DP_,GT_01_patients_with_DP_,GT_11_patients_with_DP_,GT_.._patients_with_DP_,BAF_patients_with_DP_,BAF_patients_incl_.._with_DP_)
      
      
      names(ROW)[(ncol(ROW)-6):ncol(ROW)] <- gsub("DP_",paste0("DP_",DIRECTION,"_",DP1,"_",DP2),names(ROW)[(ncol(ROW)-6):ncol(ROW)])
      
      return(ROW)
    }
    
    ROW <- BAF_BY_COVERAGE(5,5,"<=")
    ROW <- BAF_BY_COVERAGE(10,10,"<=")
    ROW <- BAF_BY_COVERAGE(11,39,"<=")
    ROW <- BAF_BY_COVERAGE(40,40,">=")
    
    
    
    
    #rearrange
    #  ROW <- ROW[,c(1:17,(ncol(ROW)-17):(ncol(ROW)),18:(ncol(ROW)-18))]
    
    ALL_ROWS <- bind_rows(ALL_ROWS,ROW)
    
  }
  
  i <- i+1
  
}
ALL_ROWS <- subset(ALL_ROWS, VQSLOD > 0)

write.table(file=OUTPUT_FILE,x=ALL_ROWS,row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')