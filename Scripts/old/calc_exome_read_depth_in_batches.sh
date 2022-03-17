#This script calculates read depths for ALL bams on mforge for a particular cancer type.

#type of cancer (e.g. LUAD, BRCA, PRAD, COAD)
CANCER_NAME=$1
MAPPED_OR_ALL=$2


#mapped: mapped reads in proper pair (-f 2), not PCR duplicate (-F 1024), with mapping quality of at least 20 (-q 20)
if [ "$MAPPED_OR_ALL" = "mapped_reads" ]; 
then SAMTOOLS_COMMAND='/home/mayo/m187735/s212975.Wickland_Immunomics/tools/samtools-1.10/bin/samtools view -c -f 2 -F 1024 -q 20';
#all: mapped and unmapped reads in proper pair (-f 1) , not PCR duplicate 
elif [ "$MAPPED_OR_ALL" = "all_reads" ]; 
then SAMTOOLS_COMMAND='/home/mayo/m187735/s212975.Wickland_Immunomics/tools/samtools-1.10/bin/samtools view -c -f 1 -F 1024';
#unmapped: unmapped reads in proper pair (-f 5), not PCR duplicate (-F 1024)
elif [ "$MAPPED_OR_ALL" = "unmapped_reads" ]; 
then SAMTOOLS_COMMAND='/home/mayo/m187735/s212975.Wickland_Immunomics/tools/samtools-1.10/bin/samtools view -c -f 5 -F 1024';

#mapped HLA: mapped reads covering HLA class I and class II region; in proper pair (-f 2), not PCR duplicate (-F 1024), with mapping quality of at least 20 (-q 20)
#first get total region: https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh38.p13
#create HLA regions file like so: grep chr6 /research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immu               nomics/gencode.v22.annotation.gtf | cut -d ';' -f1,5 | grep exon | awk '$4 >= 28510120 && $4 < 33480577' > ~/HLA.txt
#then note which genes fall at the class III interval boundaries here: https://onlinelibrary.wiley.com/cms/asset/3965cca2-f045-49d9-bc72-b3be9629df7d/tan13429-fig-0001-m.jpg
#create bed file that encompasses only the class I and class II regions
elif [ "$MAPPED_OR_ALL" = "mapped_HLA_reads" ]; 
then SAMTOOLS_COMMAND='/home/mayo/m187735/s212975.Wickland_Immunomics/tools/samtools-1.10/bin/samtools view -c -f 2 -F 1024 -q 20 -L  /research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/processing/TCGA/beds/HLA_class_I_class_II_regions.bed';

elif [ "$MAPPED_OR_ALL" = "mapped_reads_within_nimblegen_SeqCap_EZ_Exome_v2_hg38liftover_SORTED" ]; 
then SAMTOOLS_COMMAND='/home/mayo/m187735/s212975.Wickland_Immunomics/tools/samtools-1.10/bin/samtools view -c -f 2 -F 1024 -q 20 -L /research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/processing/TCGA/processed_data/capture_kit/SeqCap_EZ_Exome_v2_hg38liftover_SORTED.bed ';


elif [ "$MAPPED_OR_ALL" = "mapped_reads_within_nimblegen_SeqCap_EZ_Exome_v3_hg19_capture_targets_hg38liftover_SORTED" ]; 
then SAMTOOLS_COMMAND='/home/mayo/m187735/s212975.Wickland_Immunomics/tools/samtools-1.10/bin/samtools view -c -f 2 -F 1024 -q 20 -L /research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/processing/TCGA/processed_data/capture_kit/SeqCap_EZ_Exome_v3_hg19_capture_targets_hg38liftover_SORTED.bed ';


elif [ "$MAPPED_OR_ALL" = "mapped_reads_within_nimblegen_hg18_nimblegen_exome_version_2_hg38liftover_SORTED" ]; 
then SAMTOOLS_COMMAND='/home/mayo/m187735/s212975.Wickland_Immunomics/tools/samtools-1.10/bin/samtools view -c -f 2 -F 1024 -q 20 -L /research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/processing/TCGA/processed_data/capture_kit/hg18_nimblegen_exome_version_2_hg38liftover_SORTED.bed ';
fi;




#create directories for scripts and for outputs
if [ -d "/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/overall_read_depths/${MAPPED_OR_ALL}/${CANCER_NAME}" ]; then rm -r /home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/overall_read_depths/${MAPPED_OR_ALL}/${CANCER_NAME};
fi;

mkdir -p /home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/overall_read_depths/${MAPPED_OR_ALL}/${CANCER_NAME}

if [ -d "/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/overall_read_depths/${MAPPED_OR_ALL}/${CANCER_NAME}" ];	then rm -r /home/mayo/m187735/s212975.Wickland_Immunomics/commands_to_qsub/TCGA/${MAPPED_OR_ALL}/${CANCER_NAME};
fi;

mkdir -p /home/mayo/m187735/s212975.Wickland_Immunomics/commands_to_qsub/TCGA/${MAPPED_OR_ALL}/${CANCER_NAME}


#list all bams on mforge for cancer type
ls /research/bsi/data/controlled_access/tcga/${CANCER_NAME}/*/*gdc_realn.bam > /home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/overall_read_depths/${MAPPED_OR_ALL}/${CANCER_NAME}/${CANCER_NAME}_bam_list.txt
ls /research/bsi/data/controlled_access/tcga/${CANCER_NAME}/*/*gdc_realn_rehead.bam >> /home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/overall_read_depths/${MAPPED_OR_ALL}/${CANCER_NAME}/${CANCER_NAME}_bam_list.txt

#ls /research/bsi/data/controlled_access/tcga/downloads/wickland/cancer_types/${CANCER_NAME}/*/*gdc_realn.bam >> /home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/overall_read_depths/${MAPPED_OR_ALL}/${CANCER_NAME}/${CANCER_NAME}_bam_list.txt
#ls /research/bsi/data/controlled_access/tcga/downloads/wickland/cancer_types/${CANCER_NAME}/*/*gdc_realn_rehead.bam >> /home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/overall_read_depths/${MAPPED_OR_ALL}/${CANCER_NAME}/${CANCER_NAME}_bam_list.txt

#create scripts for sets of 100 samples
TOTAL_NUMBER_OF_BAMS=`wc -l /home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/overall_read_depths/${MAPPED_OR_ALL}/${CANCER_NAME}/${CANCER_NAME}_bam_list.txt | cut -d ' ' -f1`
NUMBER_OF_SPLITS=`echo $((TOTAL_NUMBER_OF_BAMS/100))`
BAM_LOOP_END=`echo $((NUMBER_OF_SPLITS*100))`

COUNTER=1
COUNTER2=100

#output read count (-c), filtering out unmapped (-F4), PCR duplicates (-F1024) and mapping quality below 20 (-q 20)
until [ ${COUNTER} -eq $(($BAM_LOOP_END+101)) ];  do echo 'awk -v a='"${COUNTER}"' -v b='"${COUNTER2}"' "NR==a,NR==b" /research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/processing/TCGA/overall_read_depths/'${MAPPED_OR_ALL}'/'${CANCER_NAME}'/'${CANCER_NAME}'_bam_list.txt | while read line; do echo -e -n $line"\t">> /home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/overall_read_depths/'${MAPPED_OR_ALL}'/'${CANCER_NAME}'/read_depths_samples_'${COUNTER}'-'${COUNTER2}'.txt; '${SAMTOOLS_COMMAND}' `echo $line | awk "{print $2}"` >>/home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/overall_read_depths/'${MAPPED_OR_ALL}'/'${CANCER_NAME}'/read_depths_samples_'${COUNTER}'-'${COUNTER2}'.txt;done' > /home/mayo/m187735/s212975.Wickland_Immunomics/commands_to_qsub/TCGA/${MAPPED_OR_ALL}/${CANCER_NAME}/read_depths_samples_${COUNTER}-${COUNTER2}.qsub;COUNTER=`echo $((${COUNTER}+100))`;COUNTER2=`echo $((${COUNTER}+99))`;done;


#submit qsubs
for file in /home/mayo/m187735/s212975.Wickland_Immunomics/commands_to_qsub/TCGA/${MAPPED_OR_ALL}/${CANCER_NAME}/*qsub; do qsub -l h_vmem=10G -l h_stack=50M -e /home/mayo/m187735/s212975.Wickland_Immunomics/commands_to_qsub/TCGA/${MAPPED_OR_ALL}/${CANCER_NAME}/`basename ${file} .qsub`.err -o /home/mayo/m187735/s212975.Wickland_Immunomics/commands_to_qsub/TCGA/${MAPPED_OR_ALL}/${CANCER_NAME}/`basename ${file} .qsub`.out $file;done;

