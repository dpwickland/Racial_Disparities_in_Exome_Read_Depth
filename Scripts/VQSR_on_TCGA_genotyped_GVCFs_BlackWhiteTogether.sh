#type of cancer (e.g. LUAD, BRCA, PRAD, COAD)
CANCER_NAME=$1
TISSUE=$2
PROJECT=$3


#create directory for output
OUTPUT_DIR=/research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/processing/TCGA/processed_data/VQSR_on_TCGA/${PROJECT}/${CANCER_NAME}/${TISSUE}
#rm -r $OUTPUT_DIR #delete if already exists
mkdir -p $OUTPUT_DIR
mkdir -p /home/mayo/m187735/s212975.Wickland_Immunomics/commands_to_qsub/TCGA/VQSR_on_TCGA/${PROJECT}/${CANCER_NAME}/${TISSUE}

INPUT_DIR=/research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/processing/TCGA/processed_data/GenotypeGVCFs_on_TCGA/${PROJECT}/${CANCER_NAME}/

#set VQSR variables
REF=/research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/GRCh38.d1.vd1.fa 
HAPMAP=/research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/tools/VQSR/hapmap_3.3.hg38.vcf.gz
OMNI=/research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/tools/VQSR/1000G_omni2.5.hg38.vcf.gz
G1000=/research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/tools/VQSR/1000G_phase1.snps.high_confidence.hg38.vcf.gz
DBSNP=/research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/tools/VQSR/Homo_sapiens_assembly38.dbsnp138.vcf
MILLS=/research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/tools/VQSR/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz



	
		
		INPUT=${INPUT_DIR}/${TISSUE}_BlackWhiteTogether.vcf
		
		############################
		#SNPs############
		############################
		echo "java -Xmx10g -jar /research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/tools/GenomeAnalysisTK.jar -T VariantRecalibrator -R ${REF} -input ${INPUT} -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -mode SNP -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an InbreedingCoeff -resource:hapmap,VCF,known=false,training=true,truth=true,prior=15.0 ${HAPMAP} -resource:omni,known=false,training=true,truth=true,prior=12.0 ${OMNI} -resource:1000G,known=false,training=true,truth=false,prior=10.0 ${G1000} -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${DBSNP} -recalFile ${OUTPUT_DIR}/${TISSUE}_BlackWhiteTogether.vcf_recalibrate_SNPs -tranchesFile ${OUTPUT_DIR}/${TISSUE}_BlackWhiteTogether.vcf_tranches_SNPs -rscriptFile ${OUTPUT_DIR}/${TISSUE}_BlackWhiteTogether.vcf_VQSR_SNPs_only.plots.R --disable_auto_index_creation_and_locking_when_reading_rods --target_titv 3.2
		
		java -Xmx10g -jar /research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/tools/GenomeAnalysisTK.jar -T ApplyRecalibration -R ${REF} -input ${INPUT} -o ${OUTPUT_DIR}/${TISSUE}_BlackWhiteTogether.vcf_VQSR_SNPs_only.vcf  -nt 5 -recalFile ${OUTPUT_DIR}/${TISSUE}_BlackWhiteTogether.vcf_recalibrate_SNPs -tranchesFile ${OUTPUT_DIR}/${TISSUE}_BlackWhiteTogether.vcf_tranches_SNPs --ts_filter_level 99.5 -mode SNP" > /home/mayo/m187735/s212975.Wickland_Immunomics/commands_to_qsub/TCGA/VQSR_on_TCGA/${PROJECT}/${CANCER_NAME}/${TISSUE}/BlackWhiteTogether.qsub
		
					 	
		############################
		#INDELs############
		############################
		echo "java -Xmx10g -jar /research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/tools/GenomeAnalysisTK.jar -T VariantRecalibrator -R ${REF} -input ${OUTPUT_DIR}/${TISSUE}_BlackWhiteTogether.vcf_VQSR_SNPs_only.vcf -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -mode INDEL -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an InbreedingCoeff --maxGaussians 4 -resource:mills,known=false,training=true,truth=true,prior=12.0 ${MILLS} -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${DBSNP} -recalFile ${OUTPUT_DIR}/${TISSUE}_BlackWhiteTogether.vcf_VQSR_SNPs_recalibrate_INDELs -tranchesFile ${OUTPUT_DIR}/${TISSUE}_BlackWhiteTogether.vcf_VQSR_SNPs_tranches_INDELs --disable_auto_index_creation_and_locking_when_reading_rods		
		
		
		java -Xmx10g -jar /research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/tools/GenomeAnalysisTK.jar -T ApplyRecalibration -R ${REF} -input ${OUTPUT_DIR}/${TISSUE}_BlackWhiteTogether.vcf_VQSR_SNPs_only.vcf -o ${OUTPUT_DIR}/${TISSUE}_BlackWhiteTogether_final_VQSR.vcf -recalFile ${OUTPUT_DIR}/${TISSUE}_BlackWhiteTogether.vcf_VQSR_SNPs_recalibrate_INDELs -tranchesFile ${OUTPUT_DIR}/${TISSUE}_BlackWhiteTogether.vcf_VQSR_SNPs_tranches_INDELs --ts_filter_level 99 -mode INDEL" >> /home/mayo/m187735/s212975.Wickland_Immunomics/commands_to_qsub/TCGA/VQSR_on_TCGA/${PROJECT}/${CANCER_NAME}/${TISSUE}/BlackWhiteTogether.qsub
		
		
		#submit job
		rm /research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/commands_to_qsub/TCGA/VQSR_on_TCGA/${PROJECT}/${CANCER_NAME}/${TISSUE}/BlackWhiteTogether.err 2> /dev/null
		rm /research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s212975.Wickland_Immunomics/commands_to_qsub/TCGA/VQSR_on_TCGA/${PROJECT}/${CANCER_NAME}/${TISSUE}/BlackWhiteTogether.out 2> /dev/null
		
		qsub -l h_vmem=10G -e /home/mayo/m187735/s212975.Wickland_Immunomics/commands_to_qsub/TCGA/VQSR_on_TCGA/${PROJECT}/${CANCER_NAME}/${TISSUE}/BlackWhiteTogether.err -o /home/mayo/m187735/s212975.Wickland_Immunomics/commands_to_qsub/TCGA/VQSR_on_TCGA/${PROJECT}/${CANCER_NAME}/${TISSUE}/BlackWhiteTogether.out -q lg-mem -pe threaded 6 /home/mayo/m187735/s212975.Wickland_Immunomics/commands_to_qsub/TCGA/VQSR_on_TCGA/${PROJECT}/${CANCER_NAME}/${TISSUE}/BlackWhiteTogether.qsub



