#!/bin/bash
echo "This is a script to create the report from the imputePrepSanger pipeline"

sep="========================================================================="
intro="\nThis is a report summarizing the different steps performed.\nDate: "
intro2="\n You are using version 1.1 of the pipeline:\n it uses HRCr.1 sites list and the REF allele are matched to GRCh37 \n"
today=$(date)
stringOut="$sep$intro$today$intro2 \n$sep"
resultsScreen=$1"/resultsScreen.txt"
finalReport=$1"/FinalReport.txt"

echo -e $stringOut                                                                                           >$finalReport
echo -e "The number of SNPs and individuals present in the genotype file."                                  >>$finalReport
grep "Performing single-pass .bed write (" $resultsScreen | cut -c35-                                       >>$finalReport
echo -e $sep                                                                                                >>$finalReport

#exit

echo -e "The number of SNPs flipped when using Will Rayner strand file:"                                    >>$finalReport
# This is the number of SNPs kept after matching with Will strand file
grep "SNPs flipped" $resultsScreen | grep -n 'SNPs flipped' | grep "1:.*SNPs flipped" | cut -c11-           >>$finalReport
echo -e "The number of variants remaining after removing the ones with position difference > 10bp between strand and input files:"  >>$finalReport
grep "exclude:"      $resultsScreen | grep -n 'exclude:' | grep "1:.*exclude" | cut -c14-                   >>$finalReport
echo -e "The number of variants remaining after removing the ones present in the .multiple file:"           >>$finalReport
grep "exclude:"      $resultsScreen | grep -n 'exclude:' | grep "2:.*exclude" | cut -c14-                   >>$finalReport
echo -e "The number of variants remaining (keeping only the ones in the strand file):"                      >>$finalReport
grep "extract:"     $resultsScreen | cut -c12-                                                              >>$finalReport
echo -e $sep                                                                                                >>$finalReport

# This is the number of SNPs removed in the QC steps. 
echo -e "The number of SNPs removed in the QC steps. \n"                                                    >>$finalReport
grep "(--geno)"                                          $resultsScreen                                     >>$finalReport
grep "variants removed due to Hardy-Weinberg exact test" $resultsScreen | cut -c7-                          >>$finalReport
grep "variants removed due to minor allele threshold(s)" $resultsScreen                                     >>$finalReport

echo -e $sep                                                                                                >>$finalReport
# Matching to reference panel.
echo -e "Results of the matching to reference panel: \n"                                                    >>$finalReport
echo -e "Position Matches: \n"                                                                              >>$finalReport
grep "ID matches: HRC"                         $resultsScreen                                               >>$finalReport
grep "ID Doesn't match: HRC"                   $resultsScreen                                               >>$finalReport
grep "Total Position Matches:"                 $resultsScreen                                               >>$finalReport
grep "ID Match:"                               $resultsScreen                                               >>$finalReport
grep "Different position to: HRC"              $resultsScreen                                               >>$finalReport
grep "No Match to: HRC"                        $resultsScreen                                               >>$finalReport
grep "Skipped (X, XY, Y, MT):"                 $resultsScreen                                               >>$finalReport
grep "Total in bim file:"                      $resultsScreen                                               >>$finalReport
grep "Total processed:"                        $resultsScreen                                               >>$finalReport
echo -e "\n"                                                                                                >>$finalReport
grep "Indels (ignored in r1):"                 $resultsScreen                                               >>$finalReport
echo -e "\n"                                                                                                >>$finalReport
grep "SNPs not changed:"                       $resultsScreen                                               >>$finalReport
grep "SNPs to change ref alt:"                 $resultsScreen                                               >>$finalReport
grep "Strand ok :"                             $resultsScreen                                               >>$finalReport
grep "Total Strand ok:"                        $resultsScreen                                               >>$finalReport
echo -e "\n"                                                                                                >>$finalReport
grep "Strand to change:"                       $resultsScreen                                               >>$finalReport
grep "Total checked:"                          $resultsScreen                                               >>$finalReport
grep "Total checked Strand:"                   $resultsScreen                                               >>$finalReport
grep "Total removed for allele Frequency diff" $resultsScreen                                               >>$finalReport
grep "Palindromic SNPs with Freq > 0.4 "       $resultsScreen                                               >>$finalReport
echo -e "\n"                                                                                                >>$finalReport
grep "ID and allele mismatching:"              $resultsScreen                                               >>$finalReport
grep "Non Matching alleles:"                   $resultsScreen                                               >>$finalReport
grep "Number of duplicates removed:"           $resultsScreen                                               >>$finalReport
grep "The duplicates removed can be found in"  $resultsScreen                                               >>$finalReport
echo -e "\n"                                                                                                >>$finalReport

echo -e $sep                                                                                                >>$finalReport
echo -e "The pipeline was run successfully.\n"                                                              >>$finalReport
echo -e "The resulting vcf file is named: $1/$2_afterQC-updatedChr.vcf.gz\n"                            >>$finalReport
echo -e $sep                                                                                                >>$finalReport
##


