#!/bin/bash
echo "This is a script to create the report from the imputePrepSanger pipeline"

sep="========================================================================="
intro='\nThis is a report summarizing the different steps performed.\nDate: '
today=$(date)
stringOut="$sep$intro$today \n$sep"
finalReport=$1"FinalReport.txt"

echo -e $stringOut >$finalReport
echo -e "The number of SNPs and individuals are present in the genotype file." >>$finalReport
grep  "Performing single-pass .bed write (" $1"resultsScreen.txt" | cut -c35- >>$finalReport
echo -e $sep >>finalReport

echo -e "The number of SNPs flipped when using Will Rayner strand file and the number of SNPs remaining.\n" >>$finalReport
# This is the number of SNPs kept after matching with Will strand file
grep 'SNPs flipped' $1"resultsScreen.txt" | grep -n 'SNPs flipped' | grep "1:.*SNPs flipped" | cut -c11- >>$finalReport
grep 'extract:' $1"resultsScreen.txt" | cut -c11- >>$finalReport
echo -e $sep >>$finalReport

# This is the number of SNPs removed in the QC steps. 
echo -e "The number of SNPs removed in the QC steps. \n" >>$finalReport
grep '(--geno)' $1"resultsScreen.txt" >>$finalReport
grep 'variants removed due to Hardy-Weinberg exact test' $1"resultsScreen.txt" | cut -c7- >>$finalReport
grep 'variants removed due to minor allele threshold(s)' $1"resultsScreen.txt" >>$finalReport

echo -e $sep >>$finalReport
# Matching to reference panel.
echo -e "Results of the matching to reference panel: \n" >>$finalReport
grep "Position Matches:" $1"resultsScreen.txt" >>$finalReport
grep "ID matches: HRC" $1"resultsScreen.txt" >>$finalReport
grep "ID Doesn't match: HRC" $1"resultsScreen.txt" >>$finalReport
grep "Total Position Matches:" $1"resultsScreen.txt" >>$finalReport
grep "ID Match:" $1"resultsScreen.txt" >>$finalReport
grep "Different position to: HRC" $1"resultsScreen.txt" >>$finalReport
grep "No Match to: HRC" $1"resultsScreen.txt" >>$finalReport
grep "Skipped (X, XY, Y, MT):" $1"resultsScreen.txt" >>$finalReport
grep "Total in bim file:" $1"resultsScreen.txt" >>$finalReport
grep "Total processed" $1"resultsScreen.txt" >>$finalReport
echo -e "\n" >>$finalReport
grep 'Indels (ignored in r1):' $1"resultsScreen.txt" >>$finalReport
echo -e "\n" >>$finalReport
grep 'SNPs not changed:' $1"resultsScreen.txt" >>$finalReport
grep 'SNPs to change ref alt:' $1"resultsScreen.txt" >>$finalReport
grep 'Strand ok:' $1"resultsScreen.txt" >>$finalReport
grep 'Total Strand ok :' $1"resultsScreen.txt" >>$finalReport
echo -e "\n" >>$finalReport
grep 'Strand to change:' $1"resultsScreen.txt" >>$finalReport
grep 'Total checked:' $1"resultsScreen.txt" >>$finalReport
grep 'Total checked Strand :' $1"resultsScreen.txt" >>$finalReport
grep 'Total removed for allele Frequency diff' $1"resultsScreen.txt" >>$finalReport
grep 'Palindromic SNPs with Freq > 0.4:' $1"resultsScreen.txt" >>$finalReport
echo -e "\n" >>$finalReport
grep 'ID and allele mismatching' $1"resultsScreen.txt" >>$finalReport
grep 'Non Matching alleles:' $1"resultsScreen.txt" >>$finalReport
grep 'Duplicates removed:' $1"resultsScreen.txt" >>$finalReport
echo -e "\n" >>$finalReport



echo -e $sep >>$finalReport
echo -e "The pipeline was run successfully.\n" >>$finalReport
echo -e "The resulting vcf file is named: $2_afterQC-updatedChr_vcf.vcf.gz\n" >>$finalReport
echo -e $sep >>$finalReport
##


