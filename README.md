# imputePrepSanger

This pipeline takes plink genotype files, and adjusts the strand, the positions, the reference alleles, performs quality control steps and output a vcf file that satisfies the requirement for submission to the Sanger Imputation Service (https://imputation.sanger.ac.uk/) for imputation using the Haplotype Reference Consortim reference panel. 


## Pipeline for creation of input data.

The goal of this pipeline is to facilitate the creation of input files needed for imputation using the Sanger Imputation Service. The service allows the use of the Haplotype Reference Consortium (HRC) data as a reference panel, it includes more than 32 thousand samples, from 20 different cohorts (including the UK10K and the 1000 Genomes Project).

A description of the Sanger Imputation Service can be found at the end of this document. 

Note that for now, the present pipeline creates an input file for imputation using the HRC refenrence panel only. 

### Description of the pipeline

1.  Strand are updated using the files from Will Rayner website (http://www.well.ox.ac.uk/~wrayner/strand/). There is also the option to remove variants that had more than one high quality match to the genome. 
2.  Quality Control steps (using plink):
    1.  maximum per person missing (--mind, we suggest 0.1).
    2.  maximum per SNP missing (--geno, we suggest 0.1).
    3.  Minor allele frequency (--maf, we suggest 0.05).
    4.  Hardy-Weinberg equilibrium exact test using p-value (--hwe, we suggest 5e-8).
3.  Then a perl script, modified from Will Rayner (http://www.well.ox.ac.uk/~wrayner/tools/), creates a series of plink commands to:
    1.  Update: position, ref/alt allele assignment and strand to match HRC panel.
    2.  Remove:
        1.  Variants on chromosomes: X, XY, Y and MT.
        2.  Indels.
        3.  A/T & G/C SNPs if MAF > 0.4 (palindromic SNPs, in this situation we are less confident about the strand),
        4.  SNPs with differing alleles,
        5.  No match to reference panel,
        6.  SNPs with > 0.2 allele frequency difference to the reference,
        7.  Duplicates.
4.  Create the vcf file using plink 1.9
5.  Using BCFTOOLS (version 1.3.1) rename the chromosomes to Ensembl-style chromosome names.
6.  Using BCFTOOLS run a check to make sure the reference alleles match with the ref alleles in HRC panel.
7.  Using BCFTOOLS make sure the positions are sorted (using the index function).

Note that the thresholds mentionned above at step 3 are "hard coded" (are not yet option/parameters), and therefore should be change in the perl script directly if needed. While the thresholds at step 2 are parameters to be passed in the command line, we suggest some values but these should take into account the data (number of samples, etc.).  

### How to run the pipeline

This script was created to be used on CBRAIN (https://portal.cbrain.mcgill.ca) via a docker container (see dockerhub repo: eauforest/imputeprepsanger). For this reason it assumes that the input files are in specific folders and the command line takes a large number of arguments. 


#### Softwares needed

The following softwares are needed to run the script:

1. Plink 1.9
2. BCFTOOLS (1.3.1).
3. Perl


#### Input files

The input files needed are:

1. Genotype data in Plink .map and .ped format (build 37). 
2. The corresponding .strand file (optionally also the .multiply file) from Will Rayner website: http://www.well.ox.ac.uk/~wrayner/strand/. 
3. The HRC.r1-1.GRCh37.wgs.mac5.sites.tab file (http://www.haplotype-reference-consortium.org/site)
4. The GRCh37 reference fasta (ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz)
5. A file to update the chromosome names from plink to Ensembl https://imputation.sanger.ac.uk/www/plink2ensembl.txt

Files from points 1 and 2 need to be place in the VARDATA folder, while the remaining files need to be place in the FIXDATA folder (see command line arguments). Note that these folders can have any names. 

It is important to have the version of the .strand file corresponding to the genotype array used, and to choose the version corresponding to build 37 (to match HRC reference pannel). Note that Will Rayner's files assume that the alleles are on the TOP strand. The matching .multiply file need to be added if the option to remove variants that had more than one quality match to the genome is used. It might be a good idea to remove these variants since it might means that these probes can bind to different parts of the genomes. 


#### Command line arguments

When running the script 10 arguments are needed:

1. The name of the plink files without the extension (eg.: "foo" for the files foo.map and foo.ped)
2. The name of the .strand file withitout the extension (note that this should match also the name of the .muliply file is applicable). 
3. The maximum per person missing (--mind, we suggest 0.1).
4. The maximum per SNP missing (--geno, we suggest 0.1).
5. The minor allele frequency (--maf, we suggest 0.05).
6. The Hardy-Weinberg equilibrium exact test using p-value (--hwe, we suggest 5e-8).
7. The path and name of the VARDATA folder containing the variable data (see Input Files section above).
8. The path and name of the FIXDATA folder containing the fixed data (see Input Files section above).
9. The desired name for the output folder (this folder will be created when running the pipeline).
10. A flag indicating if variants appearing in the .multiply file should be removed. Leave blank if FALSE, or write any string for TRUE. 

An example of the command line to run the pipeline with the flag to remove variants in the .multiply file.

./imputePrep_script.sh genoPlink PsychChip_15048346_B-b37 0.10 0.10 0.05 5e-8 ../new_Imputation/sangerPipeline/ressources/vardata fixData results TRUE

An example of the command line to run the pipeline without the flag to remove variants in the .multiply file.

./imputePrep_script.sh genoPlink PsychChip_15048346_B-b37 0.10 0.10 0.05 5e-8 ../new_Imputation/sangerPipeline/ressources/vardata fixData results



#### Ouput files

Two folders containing the results are created: one with the name given at point 9 above, and the other the same + "_FullOutput". The first folder will contain the resulting .vcf files to be sent to the Sanger, a text file called FinalReport.txt, which summarise the different steps, the number of variants removed and why, and a resultsScreen.txt which contains all the information appearing at the screen when running the pipeline. 

The second folder contains all intermediate files created during the pipeline. These can be deleted, but can also help undertand some details if necessary. 


## The Sanger Imputation Service

The imputation service requires input data to be in a specific format:

1. Valid VCF format.
2. All alleles on the forward strand.
3. Coordinates on GRCh37
4. REF allele matches GRCh37
5. A single VCF file, not one file per-chromosome nor per sample.
6. Records are sorted by genomic position (chromosomal order is not important).
7. Chromosomes names should be 1, 2, 3, ..., X, Y, MT.
8. If not requesting pre-phasing, then all sites and samples should be phased with no missing data.
9. Any site or sample QC should be done before uploading the data.

The Sanger Imputation Service offers the choice of using either EAGLE2 or SHAPEIT2 for the pre-phasing, and uses PBWT for imputation. There are also four choices available for the reference panel:

| Panel         | Samples          | Sites  | Comments |
| ------------- |-----------------:| ------:| -------- |
| Haplotype Reference Consortium      | 32,470 | 39M | autosomes only; SNPs only  |
| 1000 Genomes Phase 3                | 2,504  | 82M | autosomes only; SNPs, indels, complex, SVs (structural variants, genomic rearrangements larger than 50bp accounting for around 1% of the variation among human genomes).  |
| UK10K                               | 3,781  | 23M | autosomes only; biallelic SNPs and indels |
| UK10K + 1000 Genomes Phase 3        | 6,285  | 88M | autosomes only; SNPs only  |

The HRC panel should provide more accurate imputation at lower frequencies, especially in European cohorts, since it has a larger set of samples. However, it has fewer variants than the 1000 Genomes, because:

1. HRC applied a minor allele count (MAC) filter of MAC 5 (approximatively equivalent to a MAF of 0.0077%) when combining the input cohorts.
2. The 1000 Genomes Phase 3 panel includes singletons (which might not impute well).
3. HRC excluded indels.
