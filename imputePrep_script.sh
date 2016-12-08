#!/bin/bash

echo "This is a script to create vcf files for imputation on the Sanger servers"

DATASTEM=$1
STRANDFILE=$2
MIND=$3
GENO=$4
MAF=$5
HWE=$6
VARDATA=$7
FIXDATA=$8
OUTPUT=$9

PROGNAME=$(basename $0)

function error_exit
{
#	----------------------------------------------------------------
#	Function for exit due to fatal program error
#		Accepts 1 argument:
#			string containing descriptive error message
#	----------------------------------------------------------------
	echo "${PROGNAME}: ${1:-"Unknown Error"}" 1>&2
	exit 1
}

mkdir $OUTPUT

# Create binary file
echo "== Run plink, create bimnary file =="
plink --file $VARDATA/$DATASTEM --make-bed --out $OUTPUT/$DATASTEM"_binary" | tee $OUTPUT/"resultsScreen.txt"
if [ "$?" != "0" ]; then
  error_exit "Error with PLINK, creating the binary file."
fi 
 
# Update build and strand
echo "== Run update_build.sh =="
update_build.sh $OUTPUT/$DATASTEM"_binary" $VARDATA $STRANDFILE $OUTPUT/$DATASTEM"_afterAlignment" $OUTPUT | tee -a $OUTPUT/"resultsScreen.txt"
if [ "$?" != "0" ]; then
  error_exit "Error while updating the build and strand."
fi 
 
# QC steps
echo "== Run plink for QC (step 1) =="
plink --bfile $OUTPUT/$DATASTEM"_afterAlignment" --mind $MIND --geno $GENO --maf $MAF --hwe $HWE --make-bed --out $OUTPUT/$DATASTEM"_afterQC" | tee -a $OUTPUT/"resultsScreen.txt"
if [ "$?" != "0" ]; then
  error_exit "Error with PLINK, performing the QC steps."
fi 
 
# Need to perform QC before the next command.
# Also need the .bim and (from the plink --freq command) .frq files.
echo "== Run plink for QC (step 2) =="
plink --bfile $OUTPUT/$DATASTEM"_afterAlignment" --freq --out $OUTPUT/$DATASTEM"_afterAlignment_freq" | tee -a $OUTPUT/"resultsScreen.txt"
if [ "$?" != "0" ]; then
  error_exit "Error with PLINK when getting the frequency."
fi

echo "== Run HRC-1000G-check-bim_modified.pl =="
HRC-1000G-check-bim_modified.pl -b $OUTPUT/$DATASTEM"_afterAlignment.bim" -f $OUTPUT/$DATASTEM"_afterAlignment_freq.frq" -r $FIXDATA/"HRC.r1-1.GRCh37.wgs.mac5.sites.tab" -h | tee -a $OUTPUT/"resultsScreen.txt"  
if [ "$?" != "0" ]; then
  error_exit "Error while running the perl script."
fi


## Now let's create the vcf files
chmod u+x Run-plink.sh
./Run-plink.sh
if [ "$?" != "0" ]; then
  error_exit "Error while the PLINK bash file following perl script."
fi

echo "== Run bcftools =="
bcftools annotate -Oz --rename-chrs $FIXDATA/"ucsc2ensembl.txt" $OUTPUT/$DATASTEM"_afterAlignment-updated_vcf.vcf.gz" > $OUTPUT/$DATASTEM"_afterQC-updatedChr_vcf.vcf.gz"
if [ "$?" != "0" ]; then
  error_exit "Error while renaming the chromosome."
fi

cp $FIXDATA/"human_g1k_v37.fasta" $OUTPUT/ 
bcftools norm --check-ref e -f $OUTPUT/"human_g1k_v37.fasta" $OUTPUT/$DATASTEM"_afterQC-updatedChr_vcf.vcf.gz" -o $OUTPUT/$DATASTEM"_checkRef"
if [ "$?" != "0" ]; then
  error_exit "Error while checking that the REF allele matches with GRCh37 reference."
fi

bcftools index $OUTPUT/$DATASTEM"_afterQC-updatedChr_vcf.vcf.gz"
if [ "$?" != "0" ]; then
  error_exit "Error while indexing the vcf file."
fi

echo "== Run reportRedaction =="
reportRedaction.sh $OUTPUT $DATASTEM
if [ "$?" != "0" ]; then
  error_exit "Error while write report."
fi

cp -r $OUTPUT "${OUTPUT}_FullOutput"
shopt -s extglob
cd $OUTPUT
rm -rf !("${DATASTEM}_afterQC-updatedChr_vcf.vcf.gz"|FinalReport.txt|resultsScreen.txt|"${DATASTEM}_afterQC-updatedChr_vcf.vcf.gz.csi")
cd .. 

echo "The pipeline was run successfully."
