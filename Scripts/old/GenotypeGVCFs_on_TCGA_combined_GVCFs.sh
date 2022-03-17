#type of cancer (e.g. LUAD, BRCA, PRAD, COAD)
CANCER_NAME=$1
TISSUE=$2
PROJECT=$3


#create directory for output
OUTPUT_DIR=/research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/processing/TCGA/processed_data/GenotypeGVCFs_on_TCGA/${PROJECT}/${CANCER_NAME}/${TISSUE}

if [ -d "/home/mayo/m187735/s212975.Wickland_Immunomics/commands_to_qsub/TCGA/HaplotypeCaller_on_TCGA_bams/${PROJECT}/${CANCER_NAME}/${TISSUE}/" ];
then rm -r ${OUTPUT_DIR};
fi;


mkdir -p $OUTPUT_DIR
mkdir -p /home/mayo/m187735/s212975.Wickland_Immunomics/commands_to_qsub/TCGA/GenotypeGVCFs_on_TCGA/${PROJECT}/${CANCER_NAME}/${TISSUE}

INPUT_DIR=/research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/processing/TCGA/processed_data/HaplotypeCaller_on_TCGA_bams/${PROJECT}/${CANCER_NAME}/${TISSUE}/combined/

#for GROUP in {${TISSUE}_Black,${TISSUE}_White};
#	do 
	#	for GVCF in ${INPUT_DIR}/${GROUP}_samples*g.vcf; do GVCFS_LIST="${GVCFS_LIST} --variant ${GVCF}"; done
		for GVCF in ${INPUT_DIR}/*samples*g.vcf; do GVCFS_LIST="${GVCFS_LIST} --variant ${GVCF}"; done
 
		#echo "java -Xmx20g -jar /research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/tools/GenomeAnalysisTK.jar -T GenotypeGVCFs -R /research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/GRCh38.d1.vd1.fa ${GVCFS_LIST} -o ${OUTPUT_DIR}/${GROUP}.vcf --disable_auto_index_creation_and_locking_when_reading_rods -nt 24" > /home/mayo/m187735/s212975.Wickland_Immunomics/commands_to_qsub/TCGA/GenotypeGVCFs_on_TCGA/${PROJECT}/${CANCER_NAME}/${TISSUE}/${GROUP}.qsub
		echo "java -Xmx20g -jar /research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/tools/GenomeAnalysisTK.jar -T GenotypeGVCFs -R /research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/GRCh38.d1.vd1.fa ${GVCFS_LIST} -o ${OUTPUT_DIR}/BlackWhite.vcf --disable_auto_index_creation_and_locking_when_reading_rods -nt 24" > /home/mayo/m187735/s212975.Wickland_Immunomics/commands_to_qsub/TCGA/GenotypeGVCFs_on_TCGA/${PROJECT}/${CANCER_NAME}/${TISSUE}/BlackWhite.qsub
		GVCFS_LIST=""
#done


rm /research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/commands_to_qsub/TCGA/GenotypeGVCFs_on_TCGA/${PROJECT}/${CANCER_NAME}/${TISSUE}/*err 2> /dev/null 
rm /research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/commands_to_qsub/TCGA/GenotypeGVCFs_on_TCGA/${PROJECT}/${CANCER_NAME}/${TISSUE}/*out 2> /dev/null 

for file in /research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/commands_to_qsub/TCGA/GenotypeGVCFs_on_TCGA/${PROJECT}/${CANCER_NAME}/${TISSUE}/*qsub; do qsub -l h_vmem=20G -e /home/mayo/m187735/s212975.Wickland_Immunomics/commands_to_qsub/TCGA/GenotypeGVCFs_on_TCGA/${PROJECT}/${CANCER_NAME}/${TISSUE}/`basename ${file} .qsub`.err -o /home/mayo/m187735/s212975.Wickland_Immunomics/commands_to_qsub/TCGA/GenotypeGVCFs_on_TCGA/${PROJECT}/${CANCER_NAME}/${TISSUE}/`basename ${file} .qsub`.out -q 4-day -pe threaded 24 $file;done;

#Delete intermediate GVCFs to free up space
#rm -r `dirname $INPUT_DIR`/Blood_Normal