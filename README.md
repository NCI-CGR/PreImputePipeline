# PreImpute

## Imputation-Prep : Michigan/TOPMed Imputation Server
## About
The goal of this pipeline is to prepare input files needed for Imputation using Michigan/TOPMed Imputation Server.
This pipeline takes PLINK genotype files, and adjusts the strand, the positions, the reference alleles, performs quality control steps and output a vcf file for each chromosome that satisfies the requirement for submission to Michigan/TOPMed Imputation Server using HRC/1000G/TOPMed reference panels.

TOPMed server: https://imputation.biodatacatalyst.nhlbi.nih.gov/#!

HRC/1000G server: https://imputationserver.sph.umich.edu/index.html#!

Note: there are two separate snakemake scripts for TOPMed and HRC/1000G

## Description of the pipeline
1. Build and strand are updated either by liftover or using the script and files from Will Rayner website (http://www.well.ox.ac.uk/~wrayner/strand/).
	* Note: This step updates build and is optional
 	* User has an option to input either liftover chain files or use strand files from Will Rayner website
2. Quality Control steps (using plink):
	* sample call rate (default: 0.05)
	* SNP call rate (default: 0.05)
	* Minor allele frequency: (default: 0.01)
	* Hardy-Weinberg equilibrium (default: 5e-6)
3. Perl script from Will Rayner (http://www.well.ox.ac.uk/~wrayner/tools/):
	* This script is used to check a QC'd plink .bim file against the HRC, 1000G or TOPMed reference panel in advance of imputation
 	* Update:
		* position, ref/alt allele assignment and strand to match reference panel.
	* Remove:
		* A/T & G/C SNPs if MAF > 0.4 (palindromic SNPs)
		* SNPs with differing alleles
		* No match to reference panel
		* SNPs with > 0.2 allele frequency difference to the reference
		* Duplicates
4. Create the vcf files using plink

Notes:
1. Sample data provided here is in build hg19.
2. For TOPMed imputation, if the the input genotype data is in build hg19, always perform liftover/update_build to build38, since the TOPMed reference file used in Will Rayner perl script 	(step 3) in build38

## Command line arguments

module load python3;
module load slurm

1. Required:
	* ```-p```: Full path to plink files
	* ```-d```: Full path to output directory
	* ```--req_scripts```: Full path to directory with update-build.sh, HRC-1000G-check-bim.pl, liftOver (these three scripts should be placed in same directory, available in /tools/ folder)
	* ```-ref_file```: Full path to the HRC.r1-1.GRCh37.wgs.mac5.sites.tab file or 1000G 1000GP_Phase3_combined.legend.gz

	HRC.r1-1.GRCh37.wgs.mac5.sites.tab, 1000GP_Phase3_combined.legend.gz and PASS.Variants.TOPMed_freeze5_hg38_dbSNP.tab.gz files can be downloaded from https://www.chg.ox.ac.uk/~wrayner/tools/

2. Optional:
	* ```--strand_chain_file```: Full path to corresponding .strand or liftOver .chain files
		* Few liftover and strand files are provided under Ref_files/liftOver_files and Ref_files/strand_files
		* Additional strand files can be downloaded from https://www.chg.ox.ac.uk/~wrayner/strand/sourceStrand/index.html
		* liftOver files can be downloaded from
	* ```--mind```: Call rate threshold for samples
	* ```--geno```: Call rate threshold for SNPs
	* ```--maf```: Minor allele frequency threshold for SNPs
	* ```--hwe```: Hardy weinburg threshold for SNPs
	* ```--pop```: Population flag while using 1000G reference file
	* ```-u```: If pipeline was killed unexpectedly you may need this flag to rerun

3. Example usage with 1000G (no liftover or update-build)

Note: Please provide full/complete paths

```bash
cd PreImpute
python HRC_1000G/preImputeHRC.py \
	-p full_path/to/test_data_sample/input/subjects \
	-d full_path/test_data_sample/output_1000G \
	--req_scripts full_path/tools \
	--ref_file full_path/to/1000GP_Phase3_combined.legend
```

4. Example usage with TOPMed (liftover)

Note: For TOPMed imputation, if the the input genotype data is in build hg19, always perform liftover/update_build to build38, since the TOPMed reference file used is in build38

```bash
cd PreImpute
python TOPMed/preImputeTopMed.py \
	-p full_path/to/test_data_sample/input/subjects \
	-d full_path/to/test_data_sample/output_TOPMed \
	--req_scripts full_path/tools \
	--ref_file full_path/to/PASS.Variants.TOPMed_freeze5_hg38_dbSNP.tab.gz
	--strand_chain_file full_path/Ref_files/liftOver_files/hg19ToHg38.over.chain
```


## Output Folders/Files

1. build37: results 0f liftOver or  update-build.sh script
2. Plink_sub_cr – sample call rate filter files
3. Plink_snp_cr - SNP call rate filter files and updated plink files/text files after running HRC-1000G-check-bim.pl script 		(these files will have -updated extension added to file name)
4. Plink_afterQC_freq – plink frequency files after quality control steps
5.	vcf_MIS – VCF files for each chromosome ready to submit to MIS
