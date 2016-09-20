# imputePrepSanger

This pipeline takes genotype files, and adjusts the strand, the positions, the reference alleles, performs quality control steps and output a vcf file that satisfies the requirement for submittion to the Sanger Imputation Service (https://imputation.sanger.ac.uk/) for imputation using the Haplotype Reference Consortim reference panel. 


##Pipeline for creation of input data.

The goal of this pipeline is to facilitate the creation of input files needed for imputation using the Sanger Imputation Service . The service allows the use of the Haplotype Reference Consortium (HRC) data as a reference panel, it includes more than 32 thousand samples, from 20 different cohorts (including the UK10K and the 1000 Genomes Project).

A description of the Sanger Imputation Service can be found at the end of this document. 

Note that for now, the present pipeline creates input file for imputation using the HRC refenrence panel only. 

###Description of the pipeline

1.  Build and strand are updated using the files from Will Rayner website (http://www.well.ox.ac.uk/~wrayner/strand/).
2.  Quality Control steps (using plink):
    1.  maximum per person missing: 0.1 (mind).
    2.  maximum per SNP missing: 0.1 (geno).
    3.  Minor allele frequency: 0.05 (maf).
    4.  Hardy-Weinberg equilibrium exact test using p-value 5e-8 (--hwe).
3.  Then a perl script modified from Will Rayner (http://www.well.ox.ac.uk/~wrayner/tools/) creates a series of plink commands to:
    1.  Update: position, rsID, ref/alt allele assignment and strand to match HRC panel.
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

Note that all the thresholds mentionned above are "hard coded" (are not yet option/parameters), and therefore should be change in the shell script (and/or perl script) directly if needed. 

###How to run the pipeline

This script was created to use in a docker container (see docker file:...). For this reason:

* it assumes a certain folder hierarchy
* it downloads two large reference dataset to :
    * update the position, rsID, ref/alt allele assignment and strand to match HRC panel.
    * to check at the end that the REF allele in the output VCF files matches with GRCh37 reference fasta. 

If you wish to run the pipeline using the shell script, it is important to respect the folder hierarchy (or to change the path variables in the script). You might also want to download the reference dataset only once, and then comment/erase the lines in the script that download these dataset. 

####Folder Hierarchy

```
+-- data
+-- imputePrepSanger
|  +-- ressources
|    +-- strand
|    +-- HRC_refSites
|  +-- tools
|    +-- plink
|    +-- bcftools-1.3.1
|  +-- results
```

If you wish to use the shell script directly, this is where the files from github repository should be located:

1. The script imputePrep_script.sh should be located in the folder imputePrepSanger. 
2. The file "HRC-1000G-check-bim_modified.pl" should be moved to the folder imputePrepSanger/ressources/HRC_refSites/
3. The file "ucsc2ensembl.txt" should be moved to the folder imputePrepSanger/ressources/HRC_refSites/
4. The file "update_build.sh" should be moved to the folder imputePrepSanger/ressources/strand/


####Software needed

The following softwares are needed to run the script:

1. Plink 1.9
2. BCFTOOLS (1.3.1).

The script assumes that the executable, named "plink" and "bcftools", of these softwares are located in the folder "imputePrepSanger/tools/plink/" and "imputePrepSanger/tools/bcftools-1.3.1/bin/"

####Input files

The input files needed are the genotype data in Plink map and ped format and the corresponding strand file from Will Rayner website: http://www.well.ox.ac.uk/~wrayner/strand/. It is important to have the version corresponding to the genotype array used, and to choose the version corresponding to build 37. 

If running the pipeline using the shell script these input files should be in the folder data. 

When running the script two arguments are needed:

1. The name of the plink files without the extension (eg.: "foo" for the files foo.map and foo.ped)
2. The name of the strand file with the extension


##The Sanger Imputation Service

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

