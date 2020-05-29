library(ggpubr)

TMB <- read.table('/home/mayo/m187735/s212975.Wickland_Immunomics/TCGA_metadata/TCGA_mutation_load_Hoadley.txt',header=FALSE,skip=1,col.names = c('Cancer','submitter_id','submitter_id_long','silent_per_Mb','nonsilent_per_Mb'))

clinical_data_all <- data.frame()
for (CANCER_NAME in c("BRCA","UCEC","LUAD","PRAD","COAD","ACC","BLCA","SKCM","LUSC","LGG","CESC","CHOL","ESCA","GBM","HNSC","KICH","KIRP","LIHC","MESO","OV","PAAD","PCPG","READ","SARC","STAD","THYM","UCS","TGCT","UVM","LAML","THCA","KIRC")){

  #also get total number of samples per cancer (note, however, that TCGA repository claims more cases for the following: BRCA (1097 vs 1098), GBM (599 vs 617), COAD (459 vs 461), LGG (515 vs 516), LUAD (522 vs 585), OV (587 vs 608), READ (171 vs. 172), TGCT (134 vs. 150), UCEC (548 vs. 560)
  clinical_data <- read.table(paste('/home/mayo/m187735/s212975.Wickland_Immunomics/TCGA_metadata/clinical/',CANCER_NAME,'_clinical_data.txt',sep=''),header=TRUE,stringsAsFactors = FALSE)
  
  clinical_data <- clinical_data[,c("submitter_id","race")]
  clinical_data_all <- rbind(clinical_data_all,clinical_data)
  
}
rm(clinical_data)

clinical_data_all$race <- gsub("black or african american","Black",clinical_data_all$race)
clinical_data_all$race <- gsub("white","White",clinical_data_all$race)
#clinical_data_all$race <- gsub("asian","Asian",clinical_data_all$race)
clinical_data_all$race <- gsub("not reported","Unknown",clinical_data_all$race)

all_data <- merge(clinical_data_all, TMB, by='submitter_id')

all_data <- subset(all_data, (race == 'Black' | race == 'White'))
all_data <- subset(all_data, (Cancer == 'BRCA') | Cancer =='LUAD' | Cancer == 'COAD' | Cancer == 'PRAD' | Cancer == 'UCEC' | Cancer == 'KIRC' | Cancer == 'GBM' | Cancer == 'KIRP' | Cancer == 'HNSC')


####################
#PLOT EACH CANCER SEPARATELY
####################

#function for calculating n
n_fun <- function(x){return(data.frame(y=-2,label=paste0("n=",length(x),"\nmean=",round(mean(x),digits=1),"\nSE=",round(sd(x, na.rm=TRUE)/sqrt(length(x[!is.na(x)])),digits=2))))} 

all_data$race <- factor(all_data$race, levels = c('White','Black','Asian'))

ggplot(all_data, aes(x=factor(race),y=nonsilent_per_Mb,fill=race)) + xlab("") + ylab("TMB") + facet_wrap(~Cancer,nrow=2,scales='free') +
  geom_boxplot() +  labs(title=paste0("TMB in TCGA")) + 
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
  stat_summary(fun.data=n_fun,geom='text',size=3.25) + #calculate n
  coord_cartesian(ylim=c(-4,10),expand=FALSE) +

  
  stat_compare_means(size=3,method='t.test',mapping=aes(label=..p.signif..),
                     comparisons = list(c("Black","White"),c("Black","Asian"),c("White","Asian")) ,
                     label.y=c(7,8,9),#adjust position of p-value
                     tip.length = 0,
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("p<0.0001", "p<0.001", "p<0.01", "p<0.05", "NS")))


########################
#PLOT CANCERS TOGETHER
########################

all_data$race <- factor(all_data$race, levels = c('White','Black'))#,'Asian'))



############ATTEMPT 1

#ggplot(all_data, aes(x=factor(race),y=as.numeric(nonsilent_per_Mb),fill=race,color=race)) + ylab("TMB") + facet_wrap#(~Cancer,nrow=1) +
 # geom_dotplot(y=as.numeric(all_data$nonsilent_per_Mb),inherit.aes=TRUE,fill=all_data$race) +
#  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
#               geom = "crossbar", width = 0.5,color='black')



###weirdly, when I change label to 2^x, the relative levels change!!
n_fun <- function(x){return(data.frame(vjust=4,y=-7,label=paste0(length(x),"\n",round(mean(x),digits=1),"\n",round(sd(x),digits=2))))} 
#n_fun <- function(x){return(data.frame(y=-7,label=paste0(length(x),"\n",round(mean(x),digits=1),"\n",round(sd(x, na.rm=TRUE)/sqrt(length(x[!is.na(x)])),digits=2))))} 


#calculate summary statistic
Means <- aggregate(nonsilent_per_Mb+.001 ~ Cancer + race, all_data, mean)
names(Means)[3] <- 'Mean'
SDs <- aggregate(nonsilent_per_Mb+.001 ~ Cancer + race, all_data, sd)
names(SDs)[3] <- 'SD'
Ns <- aggregate(nonsilent_per_Mb+.001 ~ Cancer + race, all_data, length)
names(Ns)[3] <- 'N'
#SEs <- data.frame(cbind(as.character(Means$Cancer), as.character(Means$race),SDs$SD/sqrt(Ns$N)))
#names(SEs) <- c('Cancer','race','SE')
#SEs[,3] <- as.numeric(SEs[,3])


ggplot(all_data, aes(x=factor(race),y=log(nonsilent_per_Mb+.001,base=2),color=as.factor(race))) + facet_wrap(~Cancer, nrow=1,strip.position='bottom') + ylab("Nonsilent mutation rate (log2)") + 
  geom_boxplot(fill=NA,outlier.shape=NULL,coef=0,size=0) + 
  geom_point(aes(color=all_data$race),position=position_dodge(.9)) +
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
               geom = "crossbar", width = 0.5,color='black') +
  
  theme(strip.text = element_text(size = 12,face='bold'),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=16),
        axis.title.y = element_text(size=18,face='bold',margin=margin(t=0,r=10,b=0,l=0)),
        axis.ticks= element_blank(), 
        plot.title = element_text(hjust=0.5,size=22,face='bold'),
        legend.position = 'right',legend.direction='vertical',
        legend.text=element_text(size=16),
        plot.margin = margin(t=15,r=25,b=45,l=25),
        legend.title=element_blank(),
        panel.spacing = unit(.25, "lines")) +
  
  #scale_y_continuous(labels = function(x) format(round(2^x,digits=2))) +
  
  
  #plot table outside margins
  #vjust to put text below technical edge of plot
  coord_cartesian(clip='off') + #allow table outside plot margins
  
  geom_text(data=subset(all_data,(Cancer=='BRCA')),x=-2.50, y=-13.15, label='SAMPLE SIZE',size=3.25,color='black',hjust=0) +
  geom_text(data=subset(all_data,(Cancer=='BRCA')),x=-2.50, y=-13.90, label='MEAN',size=3.25,color='black',hjust=0) +
  geom_text(data=subset(all_data,(Cancer=='BRCA')),x=-2.50, y=-14.65, label='SD',size=3.25,color='black',hjust=0) +
  
  #geom_text(data=Means, aes(label=round(Mean,base=2),digits=2),y=-10,vjust=7)) +
  #geom_hline(data=all_data,yintercept = -0.5) +
  stat_summary(fun.data=n_fun,geom='text',size=3.25)  #calculate n






###weirdly, when I change label to 2^x, the relative levels change!!
n_fun <- function(x){return(data.frame(vjust=1.55,y=-3,label=paste0(length(x),"\n",round(mean(x),digits=1),"\n",round(sd(x),digits=2))))} 
#n_fun <- function(x){return(data.frame(y=-7,label=paste0(length(x),"\n",round(mean(x),digits=1),"\n",round(sd(x, na.rm=TRUE)/sqrt(length(x[!is.na(x)])),digits=2))))} 




ggplot(all_data, aes(x=factor(race),y=nonsilent_per_Mb,color=as.factor(race))) + facet_wrap(~Cancer, nrow=1,strip.position='bottom') + ylab("Nonsilent mutation rate") + 
  geom_boxplot(fill=NA,outlier.shape=NULL,coef=0,size=0) + 
  geom_point(aes(color=all_data$race)) +
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
               geom = "crossbar", width = 0.5,color='black') +
  
  theme(strip.text = element_text(size = 12,face='bold'),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=16),
        axis.title.y = element_text(size=18,face='bold',margin=margin(t=0,r=10,b=0,l=0)),
        axis.ticks= element_blank(), 
        plot.title = element_text(hjust=0.5,size=22,face='bold'),
        legend.position = 'right',legend.direction='vertical',
        legend.text=element_text(size=16),
        plot.margin = margin(t=15,r=25,b=45,l=25),
        legend.title=element_blank(),
        panel.spacing = unit(.25, "lines")) +
  
  #scale_y_continuous(labels = function(x) format(round(2^x,digits=2))) +
  
  
  #plot table outside margins
  #vjust to put text below technical edge of plot
  #coord_cartesian: only show plot within these bounds, but use ALL data points for calculations
  coord_cartesian(clip='off',ylim=c(0,45)) + #allow table outside plot margins
  
  geom_text(data=subset(all_data,(Cancer=='BRCA')),x=-2.50, y=-6, label='SAMPLE SIZE',size=3.25,color='black',hjust=0) +
  geom_text(data=subset(all_data,(Cancer=='BRCA')),x=-2.50, y=-7.75, label='MEAN',size=3.25,color='black',hjust=0) +
  geom_text(data=subset(all_data,(Cancer=='BRCA')),x=-2.50, y=-9.35, label='SD',size=3.25,color='black',hjust=0) +
  
  #geom_text(data=Means, aes(label=round(Mean,base=2),digits=2),y=-10,vjust=7)) +
  #geom_hline(data=all_data,yintercept = -0.5) +
  stat_summary(fun.data=n_fun,geom='text',size=3.25)  #calculate n









##########USE THIS PLOT -- HAS S JITTER

#calculate summary statistic
#Medians <- aggregate(log(nonsilent_per_Mb+.001,base=2) ~ Cancer + race, all_data, median)
Medians <- aggregate(nonsilent_per_Mb ~ Cancer + race, all_data, median)
names(Medians)[3] <- 'Median'
Medians$Median <- round(Medians$Median, digits=2)
#order medians
Medians <- Medians[order(Medians$Cancer, Medians$race),]

Medians$x <- rep(c(0.85,0.50,0.15),8) #IF AZN INCLUDED
Medians$x <- rep(c(0.85,0.50),8)
Medians$y <- -11


#SDs <- aggregate(log(nonsilent_per_Mb+.001,base=2) ~ Cancer + race, all_data, sd)
SDs <- aggregate(nonsilent_per_Mb ~ Cancer + race, all_data, sd)
names(SDs)[3] <- 'SD'
SDs$SD <- round(SDs$SD, digits=2)

Ns <- aggregate(log(nonsilent_per_Mb+.001,base=2) ~ Cancer + race, all_data, length)
names(Ns)[3] <- 'N'
#SEs <- data.frame(cbind(as.character(Means$Cancer), as.character(Means$race),SDs$SD/sqrt(Ns$N)))
#names(SEs) <- c('Cancer','race','SE')
#SEs[,3] <- as.numeric(SEs[,3])





theme_set(theme_bw()) #set theme to remove gridlines

all_data_jittered <- data.frame()
for (CANCER_NAME in c("BRCA","COAD","GBM","KIRC","KIRP","LUAD","PRAD","UCEC","HNSC")){
  one_cancer <- subset(all_data, (Cancer==CANCER_NAME))
  
  all_races_one_cancer <- data.frame()
  for (RACE in c("White","Black")){#,"Asian")){
    if (RACE == 'White'){cutoffs = c(0,0.30)}
    if (RACE == 'Black'){cutoffs = c(0.35,0.65)}
    #if (RACE == 'Asian'){cutoffs = c(0.70,1)}
    
    one_race_one_cancer <- subset(one_cancer, (race==RACE))
    one_race_one_cancer <- one_race_one_cancer[order(one_race_one_cancer$nonsilent_per_Mb,decreasing = TRUE),]#sort by nonsilent
    one_race_one_cancer$jitter <- rev(seq(cutoffs[1], cutoffs[2], length.out=(nrow(one_race_one_cancer)))) #assign jitter based on total no. rows
    all_races_one_cancer <- rbind(all_races_one_cancer, one_race_one_cancer)
  }
  
  all_data_jittered <- rbind(all_data_jittered, all_races_one_cancer)
}











ggplot(all_data_jittered, aes(y=log(nonsilent_per_Mb+.001,base=2),color=as.factor(race))) + facet_wrap(~Cancer,nrow=1,strip.position='bottom') + ylab("log2 nonsilent mutation rate") + 
  geom_point(aes(x=jitter,color=all_data_jittered$race)) +
  geom_boxplot(size=0.00005,fatten=40000,fill=NA,outlier.shape=NA,coef=0,color='black',aes(x=jitter,fill=all_data_jittered$race,group=all_data_jittered$race)) +
  #outlier.shape=NA to remove points; size to remove plot; coef=0 to rm whiskers; size=0.00005 to make box almost invisible; fatten to make median visible
  
  #  stat_summary(fun.y = Median, fun.ymin = Median, fun.ymax = Median,
  #              geom = "crossbar", width = 0.5,color='black') +
  #  scale_y_log10(breaks = c(.25, 1, 10, 100, 500)) +
  
  theme(strip.text = element_text(size = 12,face='bold'),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=16),
        axis.title.y = element_text(size=18,face='bold',margin=margin(t=0,r=10,b=0,l=0)),
        axis.ticks= element_blank(), 
        plot.title = element_text(hjust=0.5,size=22,face='bold'),
        legend.position = 'right',legend.direction='vertical',
        legend.text=element_text(size=16),
        plot.margin = margin(t=15,r=25,b=45,l=25),
        legend.title=element_blank(),
        panel.spacing = unit(.25, "lines")) +
  
  #scale_y_continuous(labels = function(x) format(round(2^x,digits=2))) +
  
  
  #plot table outside margins
  #vjust to put text below technical edge of plot
  #coord_cartesian: only show plot within these bounds, but use ALL data points for calculations
  coord_cartesian(clip='off') + #allow table outside plot margins
  
  
  #have to add each of these manually b/c normally the position is based on where the x-values fall; because of the jitter, there are MANY x-values for each cancer type
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=-1.00, y=-12.45, label='SAMPLE SIZE',size=3.25,color='black',hjust=0) +
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=-1.00, y=-13.15, label='MEDIAN',size=3.25,color='black',hjust=0) +
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=-1.00, y=-13.85, label='SD',size=3.25,color='black',hjust=0)            +
  
  
  
  #N
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=0.1, y=-12.45, label=subset(Ns,(Cancer=='BRCA' & race=='White'))$N,size=3.25,color='black',hjust=0.5)  +
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=0.46, y=-12.45, label=subset(Ns,(Cancer=='BRCA' & race=='Black'))$N,size=3.25,color='black',hjust=0.5)  + 
  #geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=0.80, y=-12.45, label=subset(Ns,(Cancer=='BRCA' & race=='Asian'))$N,size=3.25,color='black',hjust=0.5) +
  
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=1.25*.73, y=-12.45, label=subset(Ns,(Cancer=='COAD' & race=='White'))$N,size=3.25,color='black',hjust=0.5)  +
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=1.60*.73, y=-12.45, label=subset(Ns,(Cancer=='COAD' & race=='Black'))$N,size=3.25,color='black',hjust=0.5)  + 
  #geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=1.95, y=-12.45, label=subset(Ns,(Cancer=='COAD' & race=='Asian'))$N,size=3.25,color='black',hjust=0.5) + 
  
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=2.4*.70, y=-12.45, label=subset(Ns,(Cancer=='GBM' & race=='White'))$N,size=3.25,color='black',hjust=0.5)  +
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=2.80*.70, y=-12.45, label=subset(Ns,(Cancer=='GBM' & race=='Black'))$N,size=3.25,color='black',hjust=0.5)  + 
  #geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=3.1, y=-12.45, label=subset(Ns,(Cancer=='GBM' & race=='Asian'))$N,size=3.25,color='black',hjust=0.5) + 
  
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=3.55*.69, y=-12.45, label=subset(Ns,(Cancer=='HNSC' & race=='White'))$N,size=3.25,color='black',hjust=0.5)  +
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=4.0*.69, y=-12.45, label=subset(Ns,(Cancer=='HNSC' & race=='Black'))$N,size=3.25,color='black',hjust=0.5)  + 
  #geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=4.25, y=-12.45, label=subset(Ns,(Cancer=='KIRC' & race=='Asian'))$N,size=3.25,color='black',hjust=0.5) +
  
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=4.7*.69, y=-12.45, label=subset(Ns,(Cancer=='KIRC' & race=='White'))$N,size=3.25,color='black',hjust=0.5)  +
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=5.10*.69, y=-12.45, label=subset(Ns,(Cancer=='KIRC' & race=='Black'))$N,size=3.25,color='black',hjust=0.5)  + 
  #geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=5.4, y=-12.45, label=subset(Ns,(Cancer=='KIRP' & race=='Asian'))$N,size=3.25,color='black',hjust=0.5) +
  
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=5.85*.68, y=-12.45, label=subset(Ns,(Cancer=='KIRP' & race=='White'))$N,size=3.25,color='black',hjust=0.5)  +
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=6.35*.68, y=-12.45, label=subset(Ns,(Cancer=='KIRP' & race=='Black'))$N,size=3.25,color='black',hjust=0.5)  + 
  #geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=6.55, y=-12.45, label=subset(Ns,(Cancer=='LUAD' & race=='Asian'))$N,size=3.25,color='black',hjust=0.5) +
  
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=7.0*.68, y=-12.45, label=subset(Ns,(Cancer=='LUAD' & race=='White'))$N,size=3.25,color='black',hjust=0.5)  +
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=7.40*.68, y=-12.45, label=subset(Ns,(Cancer=='LUAD' & race=='Black'))$N,size=3.25,color='black',hjust=0.5)  + 
  #geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=7.7, y=-12.45, label=subset(Ns,(Cancer=='PRAD' & race=='Asian'))$N,size=3.25,color='black',hjust=0.5) +
  
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=8.15*.68, y=-12.45, label=subset(Ns,(Cancer=='PRAD' & race=='White'))$N,size=3.25,color='black',hjust=0.5)  +
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=8.60*.68, y=-12.45, label=subset(Ns,(Cancer=='PRAD' & race=='Black'))$N,size=3.25,color='black',hjust=0.5)  + 
  #geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=8.85, y=-12.45, label=subset(Ns,(Cancer=='UCEC' & race=='Asian'))$N,size=3.25,color='black',hjust=0.5) +
  
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=9.20*.68, y=-12.45, label=subset(Ns,(Cancer=='UCEC' & race=='White'))$N,size=3.25,color='black',hjust=0.5)  +
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=9.65*.68, y=-12.45, label=subset(Ns,(Cancer=='UCEC' & race=='Black'))$N,size=3.25,color='black',hjust=0.5)  + 
  #geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=8.85, y=-12.45, label=subset(Ns,(Cancer=='UCEC' & race=='Asian'))$N,size=3.25,color='black',hjust=0.5) +
  
  
  #MEDIAN
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=0.1, y=-13.15, label=subset(Medians,(Cancer=='BRCA' & race=='White'))$Median,size=3.25,color='black',hjust=0.5)  +
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=0.46, y=-13.15, label=subset(Medians,(Cancer=='BRCA' & race=='Black'))$Median,size=3.25,color='black',hjust=0.5)  + 
  #geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=0.80, y=-13.15, label=subset(Medians,(Cancer=='BRCA' & race=='Asian'))$Median,size=3.25,color='black',hjust=0.5) +
  
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=1.25*.73, y=-13.15, label=subset(Medians,(Cancer=='COAD' & race=='White'))$Median,size=3.25,color='black',hjust=0.5)  +
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=1.60*.73, y=-13.15, label=subset(Medians,(Cancer=='COAD' & race=='Black'))$Median,size=3.25,color='black',hjust=0.5)  + 
  #geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=1.95, y=-13.15, label=subset(Medians,(Cancer=='COAD' & race=='Asian'))$Median,size=3.25,color='black',hjust=0.5) + 
  
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=2.4*.70, y=-13.15, label=subset(Medians,(Cancer=='GBM' & race=='White'))$Median,size=3.25,color='black',hjust=0.5)  +
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=2.80*.70, y=-13.15, label=subset(Medians,(Cancer=='GBM' & race=='Black'))$Median,size=3.25,color='black',hjust=0.5)  + 
  #geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=3.1, y=-13.15, label=subset(Medians,(Cancer=='GBM' & race=='Asian'))$Median,size=3.25,color='black',hjust=0.5) + 
  
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=3.55*.69, y=-13.15, label=subset(Medians,(Cancer=='HNSC' & race=='White'))$Median,size=3.25,color='black',hjust=0.5)  +
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=4.0*.69, y=-13.15, label=subset(Medians,(Cancer=='HNSC' & race=='Black'))$Median,size=3.25,color='black',hjust=0.5)  + 
  #geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=4.25, y=-13.15, label=subset(Medians,(Cancer=='KIRC' & race=='Asian'))$Median,size=3.25,color='black',hjust=0.5) +
  
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=4.7*.69, y=-13.15, label=subset(Medians,(Cancer=='KIRC' & race=='White'))$Median,size=3.25,color='black',hjust=0.5)  +
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=5.10*.69, y=-13.15, label=subset(Medians,(Cancer=='KIRC' & race=='Black'))$Median,size=3.25,color='black',hjust=0.5)  + 
  #geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=5.4, y=-13.15, label=subset(Medians,(Cancer=='KIRP' & race=='Asian'))$Median,size=3.25,color='black',hjust=0.5) +
  
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=5.85*.68, y=-13.15, label=subset(Medians,(Cancer=='KIRP' & race=='White'))$Median,size=3.25,color='black',hjust=0.5)  +
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=6.35*.68, y=-13.15, label=subset(Medians,(Cancer=='KIRP' & race=='Black'))$Median,size=3.25,color='black',hjust=0.5)  + 
  #geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=6.55, y=-13.15, label=subset(Medians,(Cancer=='LUAD' & race=='Asian'))$Median,size=3.25,color='black',hjust=0.5) +
  
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=7.0*.68, y=-13.15, label=subset(Medians,(Cancer=='LUAD' & race=='White'))$Median,size=3.25,color='black',hjust=0.5)  +
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=7.40*.68, y=-13.15, label=subset(Medians,(Cancer=='LUAD' & race=='Black'))$Median,size=3.25,color='black',hjust=0.5)  + 
  #geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=7.7, y=-13.15, label=subset(Medians,(Cancer=='PRAD' & race=='Asian'))$Median,size=3.25,color='black',hjust=0.5) +
  
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=8.15*.68, y=-13.15, label=subset(Medians,(Cancer=='PRAD' & race=='White'))$Median,size=3.25,color='black',hjust=0.5)  +
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=8.60*.68, y=-13.15, label=subset(Medians,(Cancer=='PRAD' & race=='Black'))$Median,size=3.25,color='black',hjust=0.5)  + 
  #geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=8.85, y=-13.15, label=subset(Medians,(Cancer=='UCEC' & race=='Asian'))$Median,size=3.25,color='black',hjust=0.5) +
  
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=9.20*.68, y=-13.15, label=subset(Medians,(Cancer=='UCEC' & race=='White'))$Median,size=3.25,color='black',hjust=0.5)  +
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=9.65*.68, y=-13.15, label=subset(Medians,(Cancer=='UCEC' & race=='Black'))$Median,size=3.25,color='black',hjust=0.5)  + 
  #geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=8.85, y=-13.15, label=subset(Medians,(Cancer=='UCEC' & race=='Asian'))$Median,size=3.25,color='black',hjust=0.5) +
  
  
  
  
  #SD
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=0.1, y=-13.85, label=subset(SDs,(Cancer=='BRCA' & race=='White'))$SD,size=3.25,color='black',hjust=0.5)  +
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=0.46, y=-13.85, label=subset(SDs,(Cancer=='BRCA' & race=='Black'))$SD,size=3.25,color='black',hjust=0.5)  + 
  #geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=0.80, y=-13.85, label=subset(SDs,(Cancer=='BRCA' & race=='Asian'))$SD,size=3.25,color='black',hjust=0.5) +
  
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=1.25*.73, y=-13.85, label=subset(SDs,(Cancer=='COAD' & race=='White'))$SD,size=3.25,color='black',hjust=0.5)  +
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=1.60*.73, y=-13.85, label=subset(SDs,(Cancer=='COAD' & race=='Black'))$SD,size=3.25,color='black',hjust=0.5)  + 
  #geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=1.95, y=-13.85, label=subset(SDs,(Cancer=='COAD' & race=='Asian'))$SD,size=3.25,color='black',hjust=0.5) + 
  
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=2.4*.70, y=-13.85, label=subset(SDs,(Cancer=='GBM' & race=='White'))$SD,size=3.25,color='black',hjust=0.5)  +
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=2.80*.70, y=-13.85, label=subset(SDs,(Cancer=='GBM' & race=='Black'))$SD,size=3.25,color='black',hjust=0.5)  + 
  #geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=3.1, y=-13.85, label=subset(SDs,(Cancer=='GBM' & race=='Asian'))$SD,size=3.25,color='black',hjust=0.5) + 
  
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=3.55*.69, y=-13.85, label=subset(SDs,(Cancer=='HNSC' & race=='White'))$SD,size=3.25,color='black',hjust=0.5)  +
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=4.0*.69, y=-13.85, label=subset(SDs,(Cancer=='HNSC' & race=='Black'))$SD,size=3.25,color='black',hjust=0.5)  + 
  #geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=4.25, y=-13.85, label=subset(SDs,(Cancer=='KIRC' & race=='Asian'))$SD,size=3.25,color='black',hjust=0.5) +
  
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=4.7*.69, y=-13.85, label=subset(SDs,(Cancer=='KIRC' & race=='White'))$SD,size=3.25,color='black',hjust=0.5)  +
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=5.10*.69, y=-13.85, label=subset(SDs,(Cancer=='KIRC' & race=='Black'))$SD,size=3.25,color='black',hjust=0.5)  + 
  #geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=5.4, y=-13.85, label=subset(SDs,(Cancer=='KIRP' & race=='Asian'))$SD,size=3.25,color='black',hjust=0.5) +
  
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=5.85*.68, y=-13.85, label=subset(SDs,(Cancer=='KIRP' & race=='White'))$SD,size=3.25,color='black',hjust=0.5)  +
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=6.35*.68, y=-13.85, label=subset(SDs,(Cancer=='KIRP' & race=='Black'))$SD,size=3.25,color='black',hjust=0.5)  + 
  #geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=6.55, y=-13.85, label=subset(SDs,(Cancer=='LUAD' & race=='Asian'))$SD,size=3.25,color='black',hjust=0.5) +
  
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=7.0*.68, y=-13.85, label=subset(SDs,(Cancer=='LUAD' & race=='White'))$SD,size=3.25,color='black',hjust=0.5)  +
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=7.40*.68, y=-13.85, label=subset(SDs,(Cancer=='LUAD' & race=='Black'))$SD,size=3.25,color='black',hjust=0.5)  + 
  #geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=7.7, y=-13.85, label=subset(SDs,(Cancer=='PRAD' & race=='Asian'))$SD,size=3.25,color='black',hjust=0.5) +
  
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=8.15*.68, y=-13.85, label=subset(SDs,(Cancer=='PRAD' & race=='White'))$SD,size=3.25,color='black',hjust=0.5)  +
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=8.60*.68, y=-13.85, label=subset(SDs,(Cancer=='PRAD' & race=='Black'))$SD,size=3.25,color='black',hjust=0.5)  + 
  #geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=8.85, y=-13.85, label=subset(SDs,(Cancer=='UCEC' & race=='Asian'))$SD,size=3.25,color='black',hjust=0.5) +
  
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=9.20*.68, y=-13.85, label=subset(SDs,(Cancer=='UCEC' & race=='White'))$SD,size=3.25,color='black',hjust=0.5)  +
  geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=9.65*.68, y=-13.85, label=subset(SDs,(Cancer=='UCEC' & race=='Black'))$SD,size=3.25,color='black',hjust=0.5)  
  #geom_text(data=subset(all_data_jittered,(Cancer=='BRCA')),x=8.85, y=-13.85, label=subset(SDs,(Cancer=='UCEC' & race=='Asian'))$SD,size=3.25,color='black',hjust=0.5) +
  
  
  