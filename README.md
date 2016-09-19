# imputePrepSanger

The goal of this pipeline is to facilitate the creation of the input files needed for imputation using the Sanger Imputation Service (https://imputation.sanger.ac.uk/). Their service allows the use of the Haplotype Reference Consortium data as a reference panel. This very large reference panel includes more than 32 thousand samples, from 20 different cohorts (including the UK10K and the 1000 Genomes Project).

###The Sanger Imputation Service

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

##Pipeline for creation of input data.

###Description of the pipeline

1.  Build and strand are updated using the files from Will Rayner website (http://www.well.ox.ac.uk/~wrayner/strand/).
2.  QC step (using plink):
     i.  maximum per person missing: 0.1 (mind).
     ii.  maximum per SNP missing: 0.1 (geno).
  3.  Minor allele frequency: 0.01 (maf).
  4.  Hardy-Weinberg equilibrium exact test using p-value 5e-8 (--hwe).
3.  Then a perl script from Will Rayner (http://www.well.ox.ac.uk/~wrayner/tools/) creates a series of plink commands to:
     i.  Update: position, rsID, ref/alt allele assignment and strand to match HRC panel.
     ii.  Remove:
        a.  Variants on chromosomes: X, XY, Y and MT.
        b.  Indels.
        c.  A/T & G/C SNPs if MAF > 0.4 (palindromic SNPs, in this situation we are less confident about the strand),
        d.  SNPs with differing alleles,
        e.  No match to reference panel,
        f.  SNPs with > 0.2 allele frequency difference to the reference (can be removed/changed),
        g.  Duplicates.
4.  Create the vcf files using plink 1.9
5.  Using BCFTOOLS (version 1.3.1) rename the chromosomes to Ensembl-style chromosome names.
6.  Using BCFTOOLS run a check to make sure the reference alleles match with the ref alleles in HRC panel.
7.  Using BCFTOOLS make sure the positions are sorted.

###How to run the pipeline

This script was created to use in a docker container (see docker file:...). For this reason:

* it assumes a certain folder hierarchy
* it downloads two large reference dataset to :
    * update the position, rsID, ref/alt allele assignment and strand to match HRC panel.
    * to check at the end that the REF allele in the output VCF files matches with GRCh37 reference fasta. 

####Folder Hierarchy

```
+-- data
+-- imputePrepSanger
|  +-- ressources
|    +-- strand
|    +-- HRC_refSites
|  +-- tools
|    +-- plink
|    +-- bcftools
|  +-- results
```

The script imputePrep_script.sh should be located in the folder imputePrepSanger, as the other files in the github repository. 


####Software needed

The following softwares are needed to run the script:

1. Plink 1.9
2. BCFTOOLS (1.3.1).


####Input file

The input file

