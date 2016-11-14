#!/bin/bash
echo "This is a script to create vcf files for imputation on the Sanger servers"
echo | ls tools/
echo | ls ressources/
echo | ls tools/plink/
echo | ls ../data/
echo $1
echo $2

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


# Paths
dataPATH="../data/"
refStrandPATH="ressources/strand/"
intermedPATH="results/"
resultsPATH="../data/"
plinkPATH="tools/plink/"
bcftoolsPATH="tools/bcftools-1.3.1/"
hrc_RaynerCheckPATH="ressources/HRC_refSites/"

echo 'Hello, world.' >$resultsPATH"foo.txt" 
 
# Software exec
PLINK_EXEC=$plinkPATH"plink"
BCFTOOLS_EXEC=$bcftoolsPATH"bin/bcftools"
RAYNER_EXEC=$refStrandPATH"update_build.sh"

# Data file name
STRANDFILE=$dataPATH$2
#DATAFILE=$dataPATH"MAVAN_PsychChip"
#@arr = split(/./, $1);
DATASTEM=$1
DATAFILE=$dataPATH$1
 
#Unzip a few files
wget "ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz"
if [ "$?" != "0" ]; then
  error_exit "Error while downloading the files from sanger website"
fi

mv "HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz" $hrc_RaynerCheckPATH"HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz"
if [ "$?" != "0" ]; then
  error_exit "Error while moving the zip file"
fi

gunzip $hrc_RaynerCheckPATH"HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz"
if [ "$?" != "0" ]; then
  error_exit "Error while unzipping the file... "
fi
echo | ls $hrc_RaynerCheckPATH 
 
# Create binary file
$PLINK_EXEC --file $DATAFILE  --make-bed --out $intermedPATH$DATASTEM"_binary"
if [ "$?" != "0" ]; then
  error_exit "Error with PLINK, creating the binary file."
fi 
 
# Update build and strand
$RAYNER_EXEC $intermedPATH$DATASTEM"_binary" $STRANDFILE $intermedPATH$DATASTEM"_afterAlignment"
if [ "$?" != "0" ]; then
  error_exit "Error while updating the build and strand."
fi 
 
# QC steps,
#$PLINK_EXEC --bfile $intermedPATH$DATASTEM"_afterAlignment" --mind 0.1 --geno 0.1 --maf 0.05 --hwe 5e-8 --make-bed --out $intermedPATH$DATASTEM"_afterQC"
$PLINK_EXEC --bfile $intermedPATH$DATASTEM"_afterAlignment" --mind $3 --geno $4 --maf $5 --hwe $6 --make-bed --out $intermedPATH$DATASTEM"_afterQC"
if [ "$?" != "0" ]; then
  error_exit "Error with PLINK, performing the QC steps."
fi 
 
# Need to perform QC before the next command.
# Also need the .bim and (from the plink --freq command) .frq files.
$PLINK_EXEC --bfile $intermedPATH$DATASTEM"_afterQC" --freq --out $intermedPATH$DATASTEM"_afterAlignment_freq"
if [ "$?" != "0" ]; then
  error_exit "Error with PLINK when getting the frequency."
fi

perl $hrc_RaynerCheckPATH"HRC-1000G-check-bim_modified.pl" -b $intermedPATH$DATASTEM"_afterAlignment.bim" -f $intermedPATH$DATASTEM"_afterAlignment_freq.frq" -r $hrc_RaynerCheckPATH"HRC.r1-1.GRCh37.wgs.mac5.sites.tab" -h  
if [ "$?" != "0" ]; then
  error_exit "Error while running the perl script."
fi

## Now let's create the vcf files
chmod u+x Run-plink.sh
./Run-plink.sh
if [ "$?" != "0" ]; then
  error_exit "Error while the PLINK bash file following perl script."
fi
 
#Need more data for check
rm $hrc_RaynerCheckPATH"HRC.r1-1.GRCh37.wgs.mac5.sites.tab"
wget "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz"
if [ "$?" != "0" ]; then
  error_exit "Error while downloading the GRCh37 reference fasta."
fi
mv "human_g1k_v37.fasta.gz" $hrc_RaynerCheckPATH"human_g1k_v37.fasta.gz"
if [ "$?" != "0" ]; then
  error_exit "Error while moving the reference fasta file."
fi
gunzip $hrc_RaynerCheckPATH"human_g1k_v37.fasta.gz"
if [ "$?" != "0" ]; then
  error_exit "Error while unzipping the GRCh37 reference fasta file."
fi

$BCFTOOLS_EXEC annotate -Oz --rename-chrs $hrc_RaynerCheckPATH"ucsc2ensembl.txt" $intermedPATH$DATASTEM"_afterAlignment-updated_vcf.vcf.gz" > $resultsPATH$DATASTEM"_afterAlignment-updatedChr_vcf.vcf.gz"
if [ "$?" != "0" ]; then
  error_exit "Error while renaming the chromosome."
fi

$BCFTOOLS_EXEC norm --check-ref e -f $hrc_RaynerCheckPATH"human_g1k_v37.fasta" $resultsPATH$DATASTEM"_afterAlignment-updatedChr_vcf.vcf.gz" -o $resultsPATH$DATASTEM"_checkRef"
if [ "$?" != "0" ]; then
  error_exit "Error while checking that the REF allele matches with GRCh37 reference."
fi

$BCFTOOLS_EXEC index $resultsPATH$DATASTEM"_afterAlignment-updatedChr_vcf.vcf.gz"
if [ "$?" != "0" ]; then
  error_exit "Error while indexing the vcf file."
fi

echo "The pipeline was run successfully."
