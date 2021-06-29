#type of cancer (e.g. LUAD, BRCA, PRAD, COAD)
CANCER_NAME=$1
TISSUE=$2
BED_FILE=$3
PROJECT=$4

#create directories for scripts and for outputs
mkdir -p /home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/bams_lists_hg38/${PROJECT}
mkdir -p /home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/HaplotypeCaller_on_TCGA_bams/${PROJECT}/${CANCER_NAME}/${TISSUE}/White
mkdir -p /home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/HaplotypeCaller_on_TCGA_bams/${PROJECT}/${CANCER_NAME}/${TISSUE}/Black

if [ -d "/home/mayo/m187735/s212975.Wickland_Immunomics/commands_to_qsub/TCGA/HaplotypeCaller_on_TCGA_bams/${PROJECT}/${CANCER_NAME}/${TISSUE}/" ];
then rm -r /home/mayo/m187735/s212975.Wickland_Immunomics/commands_to_qsub/TCGA/HaplotypeCaller_on_TCGA_bams/${PROJECT}/${CANCER_NAME}/${TISSUE}/;
fi;

mkdir -p /home/mayo/m187735/s212975.Wickland_Immunomics/commands_to_qsub/TCGA/HaplotypeCaller_on_TCGA_bams/${PROJECT}/${CANCER_NAME}/${TISSUE}/


#list all hg38 bams on mforge for cancer type
echo " " > /home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/bams_lists_hg38/${PROJECT}/${CANCER_NAME}_bam_list_${TISSUE}.txt #first create file with blank line at top, otherwise skips when generating commands for qsub files...

#only want samples captured by a Nimblegen kit
if [ "$TISSUE" = "Primary_Tumor" ]; 
then grep -e "Black" -e "White" /research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/processing/TCGA/processed_data/overall_read_depths_and_sample_information/${CANCER_NAME}_all_reads_master_table.txt | awk '{print $1,$10,$4,$8,$14}' FS='\t' OFS='\t' | awk '$3 != "NA" && $4 == "Primary_Tumor" && $5 ~ "imble"' FS='\t' >> /home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/bams_lists_hg38/${PROJECT}/${CANCER_NAME}_bam_list_${TISSUE}.txt;

elif [ "$TISSUE" = "Blood_Normal" ]; 
then grep -e "Black" -e "White" /research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/processing/TCGA/processed_data/overall_read_depths_and_sample_information/${CANCER_NAME}_all_reads_master_table.txt | awk '{print $1,$10,$4,$8,$14}' FS='\t' OFS='\t' | awk '$3 != "NA" && $4 == "Blood_Normal" && $5 ~ "imble"' FS='\t' >> /home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/bams_lists_hg38/${PROJECT}/${CANCER_NAME}_bam_list_${TISSUE}.txt;
fi


#create scripts for sets of 10 samples
TOTAL_NUMBER_OF_BAMS=`wc -l /home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/bams_lists_hg38/${PROJECT}/${CANCER_NAME}_bam_list_${TISSUE}.txt | cut -d ' ' -f1`
NUMBER_OF_SPLITS=`echo $((TOTAL_NUMBER_OF_BAMS/10))`
BAM_LOOP_END=`echo $((NUMBER_OF_SPLITS*10))`


COUNTER=1
COUNTER2=11
#the positional variables in awk correspond to the order of cols in the bam_list.txt generated above
until [ ${COUNTER} -eq $(($BAM_LOOP_END+11)) ];  do awk -v a=${COUNTER} -v b=${COUNTER2} "NR==a,NR==b" /home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/bams_lists_hg38/${PROJECT}/${CANCER_NAME}_bam_list_${TISSUE}.txt | while read line; do awk -v PROJECT=${PROJECT} -v CANCER_NAME=${CANCER_NAME} -v BED_FILE=${BED_FILE} '{print "java -jar /research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/tools/GenomeAnalysisTK.jar -T HaplotypeCaller -R /research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/GRCh38.d1.vd1.fa -I "$3" -ERC GVCF --intervals "BED_FILE" -o  /home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/HaplotypeCaller_on_TCGA_bams/"PROJECT"/"CANCER_NAME"/"$4"/"$2"/"$1".g.vcf"}' OFS='\t' > /home/mayo/m187735/s212975.Wickland_Immunomics/commands_to_qsub/TCGA/HaplotypeCaller_on_TCGA_bams/${PROJECT}/${CANCER_NAME}/${TISSUE}/samples_${COUNTER}-${COUNTER2}.qsub;done;COUNTER=`echo $((${COUNTER}+10))`;COUNTER2=`echo $((${COUNTER}+10))`;done;

#submit qsubs
for file in /research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/commands_to_qsub/TCGA/HaplotypeCaller_on_TCGA_bams/${PROJECT}/${CANCER_NAME}/${TISSUE}/*qsub; do qsub -l h_vmem=10G -e /home/mayo/m187735/s212975.Wickland_Immunomics/commands_to_qsub/TCGA/HaplotypeCaller_on_TCGA_bams/${PROJECT}/${CANCER_NAME}/${TISSUE}/`basename ${file} .qsub`.err -o /home/mayo/m187735/s212975.Wickland_Immunomics/commands_to_qsub/TCGA/HaplotypeCaller_on_TCGA_bams/${PROJECT}/${CANCER_NAME}/${TISSUE}/`basename ${file} .qsub`.out -pe threaded 4 -q 1-day,4-day,7-day,30-day $file;done;


