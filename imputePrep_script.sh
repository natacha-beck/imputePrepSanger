#!/bin/bash
echo "This is a script to create vcf files for imputation on the Sanger servers"
echo | ls tools/
echo | ls ressources/
echo | ls tools/plink/
echo | ls ../data/
echo $1
echo $2


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
 
 
# Create binary file
$PLINK_EXEC --file $DATAFILE  --make-bed --out $intermedPATH$DATASTEM"_binary"
 
 
# Update build and strand
$RAYNER_EXEC $intermedPATH$DATASTEM"_binary" $STRANDFILE $intermedPATH$DATASTEM"_afterAlignment"
 
 
# QC steps,
$PLINK_EXEC --bfile $intermedPATH$DATASTEM"_afterAlignment" --mind 0.1 --geno 0.1 --maf 0.05 --hwe 5e-8 --make-bed --out $intermedPATH$DATASTEM"_afterQC"
 
#Unzip a few files
wget "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz"
mv "human_g1k_v37.fasta.gz" $hrc_RaynerCheckPATH"human_g1k_v37.fasta.gz"
wget "ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1/HRC.r1.GRCh37.autosomes.mac5.sites.tab.gz" 
mv "HRC.r1.GRCh37.autosomes.mac5.sites.tab.gz" $hrc_RaynerCheckPATH"HRC.r1.GRCh37.autosomes.mac5.sites.tab.gz"
gunzip $hrc_RaynerCheckPATH"human_g1k_v37.fasta.gz" 
gunzip $hrc_RaynerCheckPATH"HRC.r1.GRCh37.autosomes.mac5.sites.tab.gz"
 
 
# Need to perform QC before the next command.
# Also need the .bim and (from the plink --freq command) .frq files.
$PLINK_EXEC --bfile $intermedPATH$DATASTEM"_afterAlignment" --freq --out $intermedPATH$DATASTEM"_afterAlignment_freq"
 
perl $hrc_RaynerCheckPATH"HRC-1000G-check-bim_modified.pl" -b $intermedPATH$DATASTEM"_afterAlignment.bim" -f $intermedPATH$DATASTEM"_afterAlignment_freq.frq" -r $hrc_RaynerCheckPATH"HRC.r1-1.GRCh37.wgs.mac5.sites.tab" -h  
 
## Now let's create the vcf files
chmod u+x Run-plink.sh
./Run-plink.sh
 
$BCFTOOLS_EXEC annotate -Oz --rename-chrs $hrc_RaynerCheckPATH"ucsc2ensembl.txt" $intermedPATH$DATASTEM"_afterAlignment-updated_vcf.vcf.gz" > $resultsPATH$DATASTEM"_afterAlignment-updatedChr_vcf.vcf.gz"
 
$BCFTOOLS_EXEC norm --check-ref e -f $hrc_RaynerCheckPATH"human_g1k_v37.fasta" $resultsPATH$DATASTEM"_afterAlignment-updatedChr_vcf.vcf.gz" -o $resultsPATH$DATASTEM"_checkRef"
 
$BCFTOOLS_EXEC index $resultsPATH$DATASTEM"_afterAlignment-updatedChr_vcf.vcf.gz"
 
