#!/bin/bash
echo "This is a script to create the report from the imputePrepSanger pipeline"

sep="========================================================================="
intro='\nThis is a report summarizing the different steps performed.\nDate: '
today=$(date)
stringOut="$sep$intro$today \n$sep"

echo -e $stringOut >FinalReport.txt
echo -e "The number of SNPs and individuals are present in the genotype file." >>FinalReport.txt
grep  "Performing single-pass .bed write (" resultsScreen.txt | cut -c35- >>FinalReport.txt
echo -e $sep >>FinalReport.txt

echo -e "The number of SNPs flipped when using Will Rayner strand file and the number of SNPs remaining.\n" >>FinalReport.txt
# This is the number of SNPs kept after matching with Will strand file
grep 'SNPs flipped' resultsScreen.txt | grep -n 'SNPs flipped' | grep "1:.*SNPs flipped" | cut -c11- >>FinalReport.txt
grep 'extract:' resultsScreen.txt | cut -c11- >>FinalReport.txt
echo -e $sep >>FinalReport.txt

# This is the number of SNPs removed in the QC steps. 
echo -e "The number of SNPs removed in the QC steps. \n" >>FinalReport.txt
grep '(--geno)' resultsScreen.txt >>FinalReport.txt
grep 'variants removed due to Hardy-Weinberg exact test' resultsScreen.txt | cut -c7- >>FinalReport.txt
grep 'variants removed due to minor allele threshold(s)' resultsScreen.txt >>FinalReport.txt

echo -e $sep >>FinalReport.txt
# Matching to reference panel.
echo -e "Results of the matching to reference panel: \n" >>FinalReport.txt
grep "Position Matches:" resultsScreen.txt >>FinalReport.txt
grep "ID matches: HRC" resultsScreen.txt >>FinalReport.txt
grep "ID Doesn't match: HRC" resultsScreen.txt >>FinalReport.txt
grep "Total Position Matches:" resultsScreen.txt >>FinalReport.txt
grep "ID Match:" resultsScreen.txt >>FinalReport.txt
grep "Different position to: HRC" resultsScreen.txt >>FinalReport.txt
grep "No Match to: HRC" resultsScreen.txt >>FinalReport.txt
grep "Skipped (X, XY, Y, MT):" resultsScreen.txt >>FinalReport.txt
grep "Total in bim file:" resultsScreen.txt >>FinalReport.txt
grep "Total processed" resultsScreen.txt >>FinalReport.txt
echo -e "\n" >>FinalReport.txt
grep 'Indels (ignored in r1):' resultsScreen.txt >>FinalReport.txt
echo -e "\n" >>FinalReport.txt
grep 'SNPs not changed:' resultsScreen.txt >>FinalReport.txt
grep 'SNPs to change ref alt:' resultsScreen.txt >>FinalReport.txt
grep 'Strand ok:' resultsScreen.txt >>FinalReport.txt
grep 'Total Strand ok :' resultsScreen.txt >>FinalReport.txt
echo -e "\n" >>FinalReport.txt
grep 'Strand to change:' resultsScreen.txt >>FinalReport.txt
grep 'Total checked:' resultsScreen.txt >>FinalReport.txt
grep 'Total checked Strand :' resultsScreen.txt >>FinalReport.txt
grep 'Total removed for allele Frequency diff' resultsScreen.txt >>FinalReport.txt
grep 'Palindromic SNPs with Freq > 0.4:' resultsScreen.txt >>FinalReport.txt
echo -e "\n" >>FinalReport.txt
grep 'ID and allele mismatching' resultsScreen.txt >>FinalReport.txt
grep 'Non Matching alleles:' resultsScreen.txt >>FinalReport.txt
grep 'Duplicates removed:' resultsScreen.txt >>FinalReport.txt
echo -e "\n" >>FinalReport.txt



echo -e $sep >>FinalReport.txt
echo -e "The pipeline was run successfully.\n" >>FinalReport.txt
echo -e "The resulting vcf file is named: $1_afterQC-updatedChr_vcf.vcf.gz\n" >>FinalReport.txt
echo -e $sep >>FinalReport.txt
##

#grep  "update-chr:" resultsScreen.txt | grep -n "update-chr:" | grep  "1:.*update-chr:" >>FinalReport.txt

