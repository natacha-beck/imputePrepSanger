#!/bin/sh

#A script for updating a binary ped file using one of Will's strand files
#NRR 17th Jan 2012

#V2 13th Feb 2012. Added code to retain only SNPs in the strand file

#Required parameters:
#1. The original bed stem (not including file extension suffix)
#2. The strand file to apply
#3. The new stem for output
#Result: A new bed file (etc) using the new stem

#Required parameters:
#1. The original bed stem (not including file extension suffix)
#2. The strand file to apply
#3. The new stem for output
#Result: A new bed file (etc) using the new stem

#Added by M. Forest to match the plink version we have
PLINK_EXEC=`which plink`

#Unpack the parameters into labelled variables
stem=$1
input_dir=$2
strand_file=$3
outstem=$4
output=$5
echo Input stem is $stem
echo Strand file is $strand_file
echo Output stem is $outstem

#Cut the strand file into a file for flipping strand when necessary
flip_file=$output/$strand_file.flip
cat $input_dir/$strand_file | awk '{if ($5=="-") print $0}' | cut -f 1 > $flip_file

#Create a file that compare Will's position to position in Plink file, keep variant with difference in position less than 10bp
pos_file=$output/$strand_file.pos
awk -f checkPositions.awk $stem".bim" $input_dir/$strand_file > $pos_file

#Because Plink only allows you to update one attribute at a time, we need lots of temp
#Plink files
temp_prefix=TEMP_FILE_XX72262628_
temp1=$temp_prefix"1"

#1. Apply the flip
$PLINK_EXEC  --allow-no-sex --bfile $stem --flip $flip_file --make-bed --out $temp1
#2. Extract the SNPs in the pos file, we don't want SNPs that aren't in the strand file
$PLINK_EXEC  --allow-no-sex --bfile $temp1 --extract $pos_file --make-bed --out $outstem

#Now delete any temporary artefacts produced
rm -f $temp_prefix*

