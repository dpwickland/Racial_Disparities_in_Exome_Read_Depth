#identify shared, commonly mutated genes from 8 cancer types in cbioportal

#load data; do not read in enter line b/c contains spaces
BRCA <- read.table('/research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/cbiportal_top_genes/Mutated_Genes_BRCA.txt',skip=1)
LUAD <- read.table('/research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/cbiportal_top_genes/Mutated_Genes_LUAD.txt',skip=1)
UCEC <- read.table('/research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/cbiportal_top_genes/Mutated_Genes_UCEC.txt',skip=1)
KIRC <- read.table('/research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/cbiportal_top_genes/Mutated_Genes_KIRC.txt',skip=1)
PRAD <- read.table('/research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/cbiportal_top_genes/Mutated_Genes_PRAD.txt',skip=1)
COAD <- read.table('/research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/cbiportal_top_genes/Mutated_Genes_COAD.txt',skip=1)
GBM <- read.table('/research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/cbiportal_top_genes/Mutated_Genes_GBM.txt',skip=1)
KIRP <- read.table('/research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/cbiportal_top_genes/Mutated_Genes_KIRP.txt',skip=1)

colnames(BRCA) <- c('Gene','MutSig(Q-Value)','# mutations','Frequency','Is_cancer_gene_(source_OncoKB)')
colnames(LUAD) <- c('Gene','MutSig(Q-Value)','# mutations','Frequency','Is_cancer_gene_(source_OncoKB)')
colnames(UCEC) <- c('Gene','MutSig(Q-Value)','# mutations','Frequency','Is_cancer_gene_(source_OncoKB)')
colnames(KIRC) <- c('Gene','MutSig(Q-Value)','# mutations','Frequency','Is_cancer_gene_(source_OncoKB)')
colnames(PRAD) <- c('Gene','MutSig(Q-Value)','# mutations','Frequency','Is_cancer_gene_(source_OncoKB)')
colnames(COAD) <- c('Gene','MutSig(Q-Value)','# mutations','Frequency','Is_cancer_gene_(source_OncoKB)')
colnames(GBM) <- c('Gene','MutSig(Q-Value)','# mutations','Frequency','Is_cancer_gene_(source_OncoKB)')
colnames(KIRP) <- c('Gene','MutSig(Q-Value)','# mutations','Frequency','Is_cancer_gene_(source_OncoKB)')


#look only at gene ID, frequency, and cancer
BRCA <- BRCA[,c(1,4)]
BRCA$Cancer <- "BRCA_freq"

LUAD <- LUAD[,c(1,4)]
LUAD$Cancer <- "LUAD_freq"

UCEC <- UCEC[,c(1,4)]
UCEC$Cancer <- "UCEC_freq"

COAD <- COAD[,c(1,4)]
COAD$Cancer <- "COAD_freq"

PRAD <- PRAD[,c(1,4)]
PRAD$Cancer <- "PRAD_freq"

KIRC <- KIRC[,c(1,4)]
KIRC$Cancer <- "KIRC_freq"

KIRP <- KIRP[,c(1,4)]
KIRP$Cancer <- "KIRP_freq"

GBM <- GBM[,c(1,4)]
GBM$Cancer <- "GBM_freq"





#convert frequency to numeric with decimal, after removing any with the < or symbol
BRCA$Frequency <- gsub("<", "",BRCA$Frequency)
BRCA$Frequency <- gsub("%", "",BRCA$Frequency)
BRCA$Frequency <- as.numeric(BRCA$Frequency)
BRCA <- subset(BRCA, (Frequency > 5))

LUAD$Frequency <- gsub("<", "",LUAD$Frequency)
LUAD$Frequency <- gsub("%", "",LUAD$Frequency)
LUAD$Frequency <- as.numeric(LUAD$Frequency)
LUAD <- subset(LUAD, (Frequency > 5))

UCEC$Frequency <- gsub("<", "",UCEC$Frequency)
UCEC$Frequency <- gsub("%", "",UCEC$Frequency)
UCEC$Frequency <- as.numeric(UCEC$Frequency)
UCEC <- subset(UCEC, (Frequency > 5))

COAD$Frequency <- gsub("<", "",COAD$Frequency)
COAD$Frequency <- gsub("%", "",COAD$Frequency)
COAD$Frequency <- as.numeric(COAD$Frequency)
COAD <- subset(COAD, (Frequency > 5))

PRAD$Frequency <- gsub("<", "",PRAD$Frequency)
PRAD$Frequency <- gsub("%", "",PRAD$Frequency)
PRAD$Frequency <- as.numeric(PRAD$Frequency)
PRAD <- subset(PRAD, (Frequency > 5))

KIRC$Frequency <- gsub("<", "",KIRC$Frequency)
KIRC$Frequency <- gsub("%", "",KIRC$Frequency)
KIRC$Frequency <- as.numeric(KIRC$Frequency)
KIRC <- subset(KIRC, (Frequency > 5))

KIRP$Frequency <- gsub("<", "",KIRP$Frequency)
KIRP$Frequency <- gsub("%", "",KIRP$Frequency)
KIRP$Frequency <- as.numeric(KIRP$Frequency)
KIRP <- subset(KIRP, (Frequency > 5))

GBM$Frequency <- gsub("<", "",GBM$Frequency)
GBM$Frequency <- gsub("%", "",GBM$Frequency)
GBM$Frequency <- as.numeric(GBM$Frequency)
GBM <- subset(GBM, (Frequency > 5))


eight_cancers <- rbind(BRCA,LUAD,UCEC,KIRC,PRAD,COAD,GBM,KIRP)


#remove duplicates, which mess up the dcast fxn; keep only the gene+cancer combo with highest read depth
eight_cancers <- unique(eight_cancers, by=c('Gene','Cancer'))
eight_cancers <- aggregate(Frequency ~ Gene+Cancer, eight_cancers, max)

#count number of cancers in which gene is mutated
eight_cancers <- add_count(eight_cancers, Gene)
names(eight_cancers)[ncol(eight_cancers)] <- 'number_of_cancers_with_mutations'

#add percentage symbol
eight_cancers$Frequency <- paste0(eight_cancers$Frequency,"%")


#create table of all cancer_gene combinations
eight_cancers_table <- dcast(eight_cancers,Gene + number_of_cancers_with_mutations ~ Cancer, value.var = 'Frequency')



       