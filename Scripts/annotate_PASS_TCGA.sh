#type of cancer (e.g. LUAD, BRCA, PRAD, COAD)
PROJECT=$1
TISSUE=$2


#create directory for output
OUTPUT_DIR=/research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/processing/TCGA/processed_data/annotate_PASS_TCGA/${PROJECT}/
#rm -r $OUTPUT_DIR #delete if already exists
mkdir -p $OUTPUT_DIR
mkdir -p /home/mayo/m187735/s212975.Wickland_Immunomics/commands_to_qsub/TCGA/annotate_PASS_TCGA/${PROJECT}

INPUT_DIR=/research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/processing/TCGA/processed_data/VQSR_on_TCGA/${PROJECT}


for GROUP in {${TISSUE}_Black,${TISSUE}_White};
	do 
		
		INPUT=${INPUT_DIR}/${GROUP}_final_VQSR_chr.vcf
		
		INPUT_PASS=${OUTPUT_DIR}/`basename ${INPUT} .vcf`_PASS.vcf
		
		echo "awk '/#/ || /PASS/' ${INPUT}> ${INPUT_PASS}     		
    
		cd /research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/tools/annovar/
		
		perl /research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/tools/annovar/table_annovar.pl ${INPUT_PASS} humandb -buildver hg38 --out ${INPUT_PASS} -remove -protocol refGene,dbnsfp30a,EUR.sites.2015_08,AFR.sites.2015_08 -operation gx,f,f,f -nastring . -vcfinput " > /home/mayo/m187735/s212975.Wickland_Immunomics/commands_to_qsub/TCGA/annotate_PASS_TCGA/${PROJECT}/${CANCER_NAME}/${GROUP}.qsub;
	
done


rm /research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/commands_to_qsub/TCGA/annotate_PASS_TCGA/${PROJECT}/${CANCER_NAME}/*err 2> /dev/null
rm /research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/commands_to_qsub/TCGA/annotate_PASS_TCGA/${PROJECT}/${CANCER_NAME}/*out 2> /dev/null

for file in /research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/commands_to_qsub/TCGA/annotate_PASS_TCGA/${PROJECT}/${CANCER_NAME}/*qsub; do qsub -l h_vmem=20G -e /home/mayo/m187735/s212975.Wickland_Immunomics/commands_to_qsub/TCGA/annotate_PASS_TCGA/${PROJECT}/${CANCER_NAME}/`basename ${file} .qsub`.err -o /home/mayo/m187735/s212975.Wickland_Immunomics/commands_to_qsub/TCGA/annotate_PASS_TCGA/${PROJECT}/${CANCER_NAME}/`basename ${file} .qsub`.out -q 1-day -pe threaded 4 $file;done;



