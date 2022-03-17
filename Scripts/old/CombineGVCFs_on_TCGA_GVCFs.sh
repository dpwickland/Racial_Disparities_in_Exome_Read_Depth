#type of cancer (e.g. LUAD, BRCA, PRAD, COAD)
CANCER_NAME=$1
TISSUE=$2
PROJECT=$3


#create directory for output
mkdir /research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/processing/TCGA/processed_data/HaplotypeCaller_on_TCGA_bams/${PROJECT}/${CANCER_NAME}/${TISSUE}/combined

#list all hg38 GVCFs on mforge for cancer type
ls /home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/HaplotypeCaller_on_TCGA_bams/${PROJECT}/${CANCER_NAME}/${TISSUE}/*/*vcf > /home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/bams_lists_hg38/${PROJECT}/${CANCER_NAME}_${TISSUE}_gvcf_list.txt

grep ${TISSUE} /home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/bams_lists_hg38/${PROJECT}/${CANCER_NAME}_${TISSUE}_gvcf_list.txt | grep "White" > /home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/bams_lists_hg38/${PROJECT}/${CANCER_NAME}_gvcf_list_${TISSUE}_White.txt

grep ${TISSUE} /home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/bams_lists_hg38/${PROJECT}/${CANCER_NAME}_${TISSUE}_gvcf_list.txt | grep "Black" > /home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/bams_lists_hg38/${PROJECT}/${CANCER_NAME}_gvcf_list_${TISSUE}_Black.txt



for GROUP in {${TISSUE}_White,${TISSUE}_Black};
 do

#create scripts for sets of 50 samples
TOTAL_NUMBER_OF_GVCFS=`wc -l /home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/bams_lists_hg38/${PROJECT}/${CANCER_NAME}_gvcf_list_${GROUP}.txt | cut -d ' ' -f1`
NUMBER_OF_SPLITS=`echo $((TOTAL_NUMBER_OF_GVCFS/50))`
GVCF_LOOP_END=`echo $((NUMBER_OF_SPLITS*50))`


COUNTER=1
COUNTER2=50

 #the positional variables in awk correspond to the order of cols in the GVCF_list.txt generated above
 #group gvcfs into sets of 50, and use list of samples/paths in that group as parameter to CombineGVCFs
 #use tr to remove newline characters in the lists, replacing all except the last one with '--variant'
 
 until [ ${COUNTER} -eq $(($GVCF_LOOP_END+51)) ]; do echo -e "java -Xmx15g -jar /research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/tools/GenomeAnalysisTK.jar -T CombineGVCFs -R /research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/GRCh38.d1.vd1.fa --variant " | tr '\n' ' ' > /home/mayo/m187735/s212975.Wickland_Immunomics/commands_to_qsub/TCGA/HaplotypeCaller_on_TCGA_bams/${PROJECT}/${CANCER_NAME}/${TISSUE}/combined_${GROUP}_${COUNTER}-${COUNTER2}.qsub; awk -v a=${COUNTER} -v b=${COUNTER2} "NR==a,NR==b" /home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/bams_lists_hg38/${PROJECT}/${CANCER_NAME}_gvcf_list_${GROUP}.txt | tr '\n' '\t' | sed '$s/\t$/ /' | sed 's/\t/ --variant /g' >> /home/mayo/m187735/s212975.Wickland_Immunomics/commands_to_qsub/TCGA/HaplotypeCaller_on_TCGA_bams/${PROJECT}/${CANCER_NAME}/${TISSUE}/combined_${GROUP}_${COUNTER}-${COUNTER2}.qsub; echo "-o /home/mayo/m187735/s212975.Wickland_Immunomics/processing/TCGA/processed_data/HaplotypeCaller_on_TCGA_bams/${PROJECT}/${CANCER_NAME}/${TISSUE}/combined/${GROUP}_samples_${COUNTER}-${COUNTER2}.g.vcf" >> /home/mayo/m187735/s212975.Wickland_Immunomics/commands_to_qsub/TCGA/HaplotypeCaller_on_TCGA_bams/${PROJECT}/${CANCER_NAME}/${TISSUE}/combined_${GROUP}_${COUNTER}-${COUNTER2}.qsub;COUNTER=`echo $((${COUNTER}+50))`;COUNTER2=`echo $((${COUNTER}+49))`;done;
done


rm /research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/commands_to_qsub/TCGA/HaplotypeCaller_on_TCGA_bams/${PROJECT}/${CANCER_NAME}/${TISSUE}/*combined*err 2> /dev/null
rm /research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/commands_to_qsub/TCGA/HaplotypeCaller_on_TCGA_bams/${PROJECT}/${CANCER_NAME}/${TISSUE}/*combined*out 2> /dev/null

for file in /research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/commands_to_qsub/TCGA/HaplotypeCaller_on_TCGA_bams/${PROJECT}/${CANCER_NAME}/${TISSUE}/combined*qsub; do qsub -l h_vmem=25G -e /home/mayo/m187735/s212975.Wickland_Immunomics/commands_to_qsub/TCGA/HaplotypeCaller_on_TCGA_bams/${PROJECT}/${CANCER_NAME}/${TISSUE}/`basename ${file} .qsub`.err -o /home/mayo/m187735/s212975.Wickland_Immunomics/commands_to_qsub/TCGA/HaplotypeCaller_on_TCGA_bams/${PROJECT}/${CANCER_NAME}/${TISSUE}/`basename ${file} .qsub`.out -q 1-day,4-day $file;done;

