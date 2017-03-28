# imputePrepSanger

This pipeline takes plink genotype files, and adjusts the strand, the positions, the reference alleles, performs quality control steps and output a vcf file that satisfies the requirement for submittion to the Sanger Imputation Service (https://imputation.sanger.ac.uk/) for imputation using the Haplotype Reference Consortim reference panel. 


## Pipeline for creation of input data.

The goal of this pipeline is to facilitate the creation of input files needed for imputation using the Sanger Imputation Service . The service allows the use of the Haplotype Reference Consortium (HRC) data as a reference panel, it includes more than 32 thousand samples, from 20 different cohorts (including the UK10K and the 1000 Genomes Project).

A description of the Sanger Imputation Service can be found at the end of this document. 

Note that for now, the present pipeline creates input file for imputation using the HRC refenrence panel only. 

### Description of the pipeline

1.  Strand are updated using the files from Will Rayner website (http://www.well.ox.ac.uk/~wrayner/strand/). There is also the option to remove variants that had more than one high quality match to the genome. 
2.  Quality Control steps (using plink):
    1.  maximum per person missing: 0.1 (mind).
    2.  maximum per SNP missing: 0.1 (geno).
    3.  Minor allele frequency: 0.05 (maf).
    4.  Hardy-Weinberg equilibrium exact test using p-value 5e-8 (--hwe).
3.  Then a perl script modified from Will Rayner (http://www.well.ox.ac.uk/~wrayner/tools/) creates a series of plink commands to:
    1.  Update: position, ref/alt allele assignment and strand to match HRC panel.
    2.  Remove:
        1.  Variants on chromosomes: X, XY, Y and MT.
        2.  Indels.
        3.  A/T & G/C SNPs if MAF > 0.4 (palindromic SNPs, in this situation we are less confident about the strand),
        4.  SNPs with differing alleles,
        5.  No match to reference panel,
        6.  SNPs with > 0.2 allele frequency difference to the reference,
        7.  Duplicates.
4.  Create the vcf files using plink 1.9
5.  Using BCFTOOLS (version 1.3.1) rename the chromosomes to Ensembl-style chromosome names.
6.  Using BCFTOOLS run a check to make sure the reference alleles match with the ref alleles in HRC panel.
7.  Using BCFTOOLS make sure the positions are sorted.

Note that the thresholds mentionned above at step 3 are "hard coded" (are not yet option/parameters), and therefore should be change in the perl script directly if needed. 

### How to run the pipeline

This script was created to use in a docker container (see docker file:...). For this reason:

* it assumes a certain folder hierarchy

If you wish to run the pipeline using the shell script, it is important to respect the folder hierarchy (or to change the path variables in the script). You might also want to download the reference dataset only once, and then comment/erase the lines in the script that download these dataset. 

#### Files needed

If you wish to use the shell script directly, this is where the files from github repository should be located:

1. The file "ucsc2ensembl.txt" should be moved to the folder imputePrepSanger/ressources/HRC_refSites/


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

Files from points 1 and 2 need to be place in the $VARDATA folder, while the remaining files need to be place in the $FIXDATA folder (see command line arguments). 

It is important to have the version of the .strand file corresponding to the genotype array used, and to choose the version corresponding to build 37 (to match HRC reference pannel). Note that Will Rayner's files assume that the alleles are on the TOP strand.. The matching .multiply file need to be added if the option to remove variants that had more than one quality match to the genome is used. 



When running the script two arguments are needed:

1. The name of the plink files without the extension (eg.: "foo" for the files foo.map and foo.ped)
2. The name of the .strand file withitout the extension (note that this should match also the name of the .muliply file is applicable). 


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
