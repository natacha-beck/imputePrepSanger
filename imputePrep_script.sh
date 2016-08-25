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
refStrandPATH="../data/"
resultsPATH="../data/"
plinkPATH="tools/plink/"
bcftoolsPATH="tools/bcftools/"
hrc_RaynerCheckPATH="ressources/HRC_refSites/"

echo 'Hello, world.' >$resultsPATH"foo.txt" 
 
# Software exec
PLINK_EXEC=$plinkPATH"plink"
BCFTOOLS_EXEC=$bcftoolsPATH"bin/bcftools"
RAYNER_EXEC=$refStrandPATH"update_build.sh"

# Data file name
STRANDFILE=$refStrandPATH$2
#DATAFILE=$dataPATH"MAVAN_PsychChip"
#@arr = split(/./, $1);
DATASTEM=$1
DATAFILE=$1
 
 
# Create binary file
$PLINK_EXEC --file $results$DATAFILE  --make-bed --out $resultsPATH$DATASTEM"_binary"
 
 
# Update build and strand
$RAYNER_EXEC $resultsPATH$DATASTEM"_binary" $STRANDFILE $resultsPATH$DATASTEM"_afterAlignment"
 
 
# QC steps,
$PLINK_EXEC --bfile $resultsPATH$DATASTEM"_afterAlignment" --mind 0.1 --geno 0.1 --maf 0.05 --hwe 5e-8 --make-bed --out $resultsPATH$DATASTEM"_afterQC"
 
#Unzip a few files
wget "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz"
mv "human_g1k_v37.fasta.gz" "ressources/HRC_refSites/human_g1k_v37.fasta.gz"
wget "ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1/HRC.r1.GRCh37.autosomes.mac5.sites.tab.gz" 
mv "HRC.r1.GRCh37.autosomes.mac5.sites.tab.gz" "ressources/HRC_refSites/HRC.r1.GRCh37.autosomes.mac5.sites.tab.gz"
gunzip "ressources/HRC_refSites/human_g1k_v37.fasta.gz" 
gunzip $hrc_RaynerCheckPATH"HRC.r1-1.GRCh37.wgs.mac5.sites.tab"
 
 
# Need to perform QC before the next command.
# Also need the .bim and (from the plink --freq command) .frq files.
$PLINK_EXEC --bfile $resultsPATH$DATASTEM"_afterAlignment" --freq --out $resultsPATH$DATASTEM"_afterAlignment_freq"
 
perl $hrc_RaynerCheckPATH"HRC-1000G-check-bim_modified.pl" -b $resultsPATH$DATASTEM"_afterAlignment.bim" -f $resultsPATH$DATASTEM"_afterAlignment_freq.frq" -r $hrc_RaynerCheckPATH"HRC.r1-1.GRCh37.wgs.mac5.sites.tab" -h  
 
## Now let's create the vcf files
chmod u+x Run-plink.sh
./Run-plink.sh
 
$BCFTOOLS_EXEC annotate -Oz --rename-chrs ressources/HRC_refSites/ucsc2ensembl.txt results/$DATASTEM"_afterAlignment-updated_vcf.vcf.gz" > results/$DATASTEM"_afterAlignment-updatedChr_vcf.vcf.gz"
 
$BCFTOOLS_EXEC norm --check-ref e -f ressources/HRC_refSites/human_g1k_v37.fasta results/$DATASTEM"_afterAlignment-updatedChr_vcf.vcf.gz" -o $DATASTEM"_checkRef"
 
$BCFTOOLS_EXEC index results/$DATASTEM"_afterAlignment-updatedChr_vcf.vcf.gz"
 
