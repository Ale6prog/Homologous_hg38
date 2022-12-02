# Chameleolyser38

✍ Authors:

**BEL Alexis**

## 🤔 Context

Chameleolyser is a bioinformatics tool to identify genetic variants in homologous regions using whole-exome sequencing (WES) data. You can find the originl git on https://github.com/Genome-Bioinformatics-RadboudUMC/Chameleolyser with work on the genome hg19.
These variants remain hidden in a regular WES analysis. This implementation of our software is hg38-based.

## 🚀 Launching this program

### 🐍 Conda environment

It is highly recommended to install all dependencies by cloning the Chameleolyser repository onto your machine. 
```
git clone https://github.com/Ale6prog/Homologous_hg38
cd Chameleolyser/
conda env create -f ChameleolyserEnvironment.yml
conda activate Chameleolyser
```
### ⚙️ Usage

### Optional Prepare BED

The prepareBED function will download all necessary BED file but if you have already clone the git. Move the BED file in your working directories.
The working directory is the directory in which all intermediate and result files will be written. Choose an existing directory for this. The prepareBED function only need to be run once (also in case multiple samples are analysed in the same working directory).
```
perl Chameleolyser.pl --PrepareBED --WORKING_DIR=<WORKING_DIRECTORY>
```
### Mask reference genome

The MaskReferenceGenome function will download a copy of the hg38 reference genome. After completion, it will create a masked version of it. This option only need to be run once (also in case multiple samples are analysed in the same working directory).
```
perl Chameleolyser.pl --MaskReferenceGenome --WORKING_DIR=<WORKING_DIRECTORY>
```

### Generate masked alignments and raw VCF

This function will extract reads in the homologous regions and re-align them to the masked reference sequence. Subsequently it will call variants with a sensitive method. The sample name is an identifier of choice. The alignment filepath is the full path of the CRAM/BAM file of your sample of interest which is stored on your machine.

```
perl Chameleolyser.pl --GenerateMaskedAlignmentAndVcf --WORKING_DIR=<WORKING_DIRECTORY> --PREFIX=chr --SAMPLE_NAME=<SAMPLE_NAME> --ALIGNMENT_FP=<ALIGNMENT_FP> --NR_OF_THREADS=<NR_OF_THREADS>
```
### Filter raw variants

```
perl Chameleolyser.pl --FilterRawVariants --WORKING_DIR=<WORKING_DIRECTORY> --PREFIX=chr --SAMPLE_NAME=<SAMPLE_NAME>
```
















