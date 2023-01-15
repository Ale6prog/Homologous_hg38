# Chameleolyser38

✍ Authors:

**BEL Alexis**

## 🤔 Context

Chameleolyser is a bioinformatics tool to identify genetic variants in homologous regions using whole-exome sequencing (WES) data. You can find the originl git on https://github.com/Genome-Bioinformatics-RadboudUMC/Chameleolyser with work on the genome hg19.
These variants remain hidden in a regular WES analysis. This implementation of our software is hg38-based.

## BED

The files from the BED folder come from https://github.com/Genome-Bioinformatics-RadboudUMC/ChameleolyserBEDs. They were liftover from hg19 to HG38 with the online tool Lift Genome Annotations (https://genome.ucsc.edu/cgi-bin/hgLiftOver). It contains also another folder named error that contains all the error from the liftover. The file error first and second originate from the files which are not bed files and that contains two positions.

### BED files

#### use_masking_hg38.bed
Contains the position of one member of a group of homologous

#### use_extraction_hg38.bed
Contains the position of all homologous regions with 200 bases more.

#### use_calling_hg38.bed
Contains the true positions of all homologous regions

#### use_homologous_exons.bed
Contains the exons of the homologous regions.

#### use_all_sd.bed and use_mappedsd.txt
use_all_sd.bed is a list of sequence differences (between all homologs used by the software). 
use_mappedsd.txt are couples of SDs. If 2 sequences A and B are nearly identical then A is different relative to B and B is different relative to A, so you have 2 differences which are coupled. XXXAXXX versus XXXBXXX then you have A>B from the perspective of the first sequence and you have B>A from the perspective of the last sequence. 

#### PosToRegionID_hg38.txt.gz and RegionID_ToStrand.txt
Contains the ID for each positions and the diretion of the gene.

#### use_Sitetoexclude_hg38.bed
In regions where you locally have a lot of SDs, there will, on average, be more false positive variants. For that reason we iterated over the SDs and when 5 SDs are in a range of 10 nucleotides we excluded this subregion.

#### use
Postion of interrest

#### use_HP_hg38.bed.gz
In regions of homopolymers we identified quite some false positives due to mapping difficulties and sequencing errors. So this is again to keep accuracy high. https://github.com/ga4gh/benchmarking-tools

#### CohortAlleleFreq_chr_hg38.txt
Contains the frequencies of variations in a populations.

## 🚀 Launching this program

### 🐍 Conda environment

It is highly recommended to install all dependencies by cloning the Chameleolyser repository onto your machine. 
```
git clone https://github.com/Ale6prog/Homologous_hg38
cd Homologous_hg38/
conda env create -f ChameleolyserEnvironment.yml
conda activate Chameleolyser
```
### ⚙️ Usage

### Optional Prepare BED

The prepareBED function will download all necessary BED file but if you have already clone the git. Move the BED file in your working directories.
The working directory is the directory in which all intermediate and result files will be written. Choose an existing directory for this. The prepareBED function only need to be run once (also in case multiple samples are analysed in the same working directory). Also make sure that all file are not zipped. The program cannot read zipped file and will output empty filtered files.
```
perl Chameleolyser.pl --PrepareBED --WORKING_DIR=<WORKING_DIRECTORY>
```
### Mask reference genome

The MaskReferenceGenome function will download a copy of the hg38 reference genome. After completion, it will create a masked version of it. This option only need to be run once (also in case multiple samples are analysed in the same working directory).
```
perl Chameleolyser.pl --MaskReferenceGenome --WORKING_DIR=<WORKING_DIRECTORY>
```

### Generate masked alignments and raw VCF

This function will extract reads in the homologous regions and re-align them to the masked reference sequence. Subsequently it will call variants with a sensitive method. The sample name is an identifier of choice. The alignment filepath is the full path of the CRAM/BAM file of your sample of interest which is stored on your machine. The BAM file must have been done with the same genome of reference as the one given by the program!

```
perl Chameleolyser.pl --GenerateMaskedAlignmentAndVcf --WORKING_DIR=<WORKING_DIRECTORY> --SAMPLE_NAME=<SAMPLE_NAME> --ALIGNMENT_FP=<ALIGNMENT_FP> --NR_OF_THREADS=<NR_OF_THREADS>
```
### Filter raw variants

```
perl Chameleolyser.pl --FilterRawVariants --WORKING_DIR=<WORKING_DIRECTORY> --SAMPLE_NAME=<SAMPLE_NAME>
```
















