####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
########## 	  	 		   						     	 	 					       	  				  ########## 
########## 	  	 		   					   		  CHAMELEOLYSER  	 	 				  	 		  ##########
########## 	  	 		   						     	 	 					       	  				  ##########
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################

##################################################################################################
##################################################################################################
##################################################################################################
#######     		     			      	Load Libraries   	  		   	 			   #######
##################################################################################################
##################################################################################################
##################################################################################################

use Getopt::Long;
use IO::Zlib;
use strict;
use warnings;
use Data::Dumper;
use List::Util qw(min max);

##################################################################################################
##################################################################################################
##################################################################################################
#######     		     		 	  	  Initialize & Declare  	  		   			   #######
##################################################################################################
##################################################################################################
##################################################################################################

######################################################
### 				PROGRAM SETTINGS		   	   ###
######################################################

my 	$WORKING_DIR									= "";
my 	$SAMPLE_NAME									= "";
my 	$ALIGNMENT_FP									= "";
my	$NR_OF_THREADS									= 1;

######################################################
### 				PROGRAM OPTIONS		   	   	   ###
######################################################

my 	$PrepareBED										= 0;
my 	$MaskReferenceGenome							= 0;
my 	$GenerateMaskedAlignmentAndVcf					= 0;
my 	$FilterRawVariants								= 0;

##################################################################################################
##################################################################################################
##################################################################################################
#######     		     		  	   	  	 Get Options  	  		   			  	       #######
##################################################################################################
##################################################################################################
##################################################################################################

GetOptions ("WORKING_DIR=s"							=>	\$WORKING_DIR,
			"SAMPLE_NAME=s"							=>	\$SAMPLE_NAME,
			"ALIGNMENT_FP=s"						=>	\$ALIGNMENT_FP,
			"NR_OF_THREADS=i"						=>	\$NR_OF_THREADS,
			"PrepareBED"							=>	\$PrepareBED,
			"MaskReferenceGenome"					=>	\$MaskReferenceGenome,
			"GenerateMaskedAlignmentAndVcf"			=>	\$GenerateMaskedAlignmentAndVcf,
			"FilterRawVariants"						=>	\$FilterRawVariants);

##################################################################################################
##################################################################################################
##################################################################################################
#######     		     		  	  	  	 MAIN PROGRAM  	  		   				   	   #######
##################################################################################################
##################################################################################################
##################################################################################################

PrepareBED									($WORKING_DIR
											 ) 							if $PrepareBED;

MaskReferenceGenome 						($WORKING_DIR,
											) 							if $MaskReferenceGenome;
											 
GenerateMaskedAlignmentAndVcf				($WORKING_DIR,
											 $SAMPLE_NAME,
											 $ALIGNMENT_FP,
											 $NR_OF_THREADS) 					if $GenerateMaskedAlignmentAndVcf;
											 
FilterRawVariants							($WORKING_DIR,
											 $SAMPLE_NAME)						if $FilterRawVariants;

##################################################################################################
##################################################################################################
##################################################################################################
#######     		     			   	   	  SUBROUTINES  	  		   				 	   #######
##################################################################################################
##################################################################################################
##################################################################################################


sub PrepareBED {
	
	(my $WORKING_DIR)
	= @_;
	
	# test if working directory exists
	
	if (! -d $WORKING_DIR){
		print "$WORKING_DIR should already exist.\n";
		print "Please enter an existing working directory.\n";
		die("$!\n");
	}
	
	# download BEDs
	
	unless	(-d "$WORKING_DIR/BED/"){system("mkdir $WORKING_DIR/BED/");}
	chdir	("$WORKING_DIR/BED/");
		
	system	("wget https://raw.githubusercontent.com/Genome-Bioinformatics-RadboudUMC/ChameleolyserBEDs/main/All.formasking.noalt.chr.bed");
	system	("wget https://raw.githubusercontent.com/Genome-Bioinformatics-RadboudUMC/ChameleolyserBEDs/main/All.forextraction.noalt.chr.bed");
	system	("wget https://raw.githubusercontent.com/Genome-Bioinformatics-RadboudUMC/ChameleolyserBEDs/main/All.forvarcall.noalt.chr.bed");
	system	("wget https://raw.githubusercontent.com/Genome-Bioinformatics-RadboudUMC/ChameleolyserBEDs/main/All.homologousexons.noalt.chr.bed");
	system	("wget https://raw.githubusercontent.com/Genome-Bioinformatics-RadboudUMC/ChameleolyserBEDs/main/AllSDs.noalt.chr.bed");
	system	("wget https://raw.githubusercontent.com/Genome-Bioinformatics-RadboudUMC/ChameleolyserBEDs/main/AllToMapOnToOther.chr.txt.gz");
	system	("wget https://raw.githubusercontent.com/Genome-Bioinformatics-RadboudUMC/ChameleolyserBEDs/main/MappedSDs.chr.txt");
	system	("wget https://raw.githubusercontent.com/Genome-Bioinformatics-RadboudUMC/ChameleolyserBEDs/main/PosToRegionID.chr.txt.gz");
	
	system	("wget https://raw.githubusercontent.com/Genome-Bioinformatics-RadboudUMC/ChameleolyserBEDs/main/HPs.noalt.chr.bed.gz");
	system	("wget https://raw.githubusercontent.com/Genome-Bioinformatics-RadboudUMC/ChameleolyserBEDs/main/SitesToExclude.noalt.chr.bed");
	system	("wget https://raw.githubusercontent.com/Genome-Bioinformatics-RadboudUMC/ChameleolyserBEDs/main/CohortAlleleFreq.chr.txt");
		
	system("gunzip HPs.noalt.chr.bed.gz");

	system	("wget https://raw.githubusercontent.com/Genome-Bioinformatics-RadboudUMC/ChameleolyserBEDs/main/RegionID_ToStrand.txt");
}

sub MaskReferenceGenome {
	
	(my $WORKING_DIR)
	= @_;
	
	# find Picard
	
	my 	$PicardJarPath = `which picard`;
		$PicardJarPath =~ s/\n//g;
		$PicardJarPath =~ s/bin\/.*$/bin\//g;
		$PicardJarPath .= "../share/picard-2.20.8-0/picard.jar";
	
	# test if working directory exists
	
	if (! -d $WORKING_DIR){
		print "$WORKING_DIR should already exist.\n";
		print "Please enter an existing working directory.\n";
		die("$!\n");
	}
	
	# make new subdirs
	
	unless(-d "$WORKING_DIR/REF/"){system("mkdir $WORKING_DIR/REF/");}
		
	# download RefSeq #
	
	chdir	("$WORKING_DIR/REF/");
	print "Starting downloading hg38 patch 14 \n";
		# download data
	system	("wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.40_GRCh38.p14/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.gz");
		# decompress data
	system 	("gunzip GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.gz");
		# Change name of the file
	system	("mv GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna hg38.chr.fa");
		# Index the genome > burrows wheeler compressed index
	system 	("bwa index hg38.chr.fa");
		#  Index the region of the genome as an .fai
		# contains lenght ,position in the file; nmes length of the fasta files
	system 	("samtools faidx hg38.chr.fa");
	die;
	print "Downloading and indexing of hg38 patch 14 finish\n";

		# Recuperate the position and the lenght in the file
	system	("awk -v OFS='\t' {'print \$1,\$2'} hg38.chr.fa.fai > hg38.chr.txt");
		# Create dictionnary with the sequencé of the genome 
	system	("java -Xmx8G -jar $PicardJarPath CreateSequenceDictionary REFERENCE=hg38.chr.fa OUTPUT=hg38.chr.dict");
		
		# masks sequences in genome based on intervals defined in a feature file (bed)
	print "Starting the masking of hg38 patch 14\n";
	system 	("bedtools maskfasta -fi hg38.chr.fa -bed ../BED/All.formasking.chr.hg38.bed -fo hg38.masked.chr.fa");
		# Index the masked genome
	system 	("bwa index hg38.masked.chr.fa");
		# Index the region of the masked genome as a fai file
	system 	("samtools faidx hg38.masked.chr.fa");	
		# Create dictionnary with the sequence of the genome 
	system	("java -Xmx8G -jar $PicardJarPath CreateSequenceDictionary REFERENCE=hg38.masked.chr.fa OUTPUT=hg38.masked.chr.dict");
	print "Masking of hg38 patch 14 Ready!\n";
}


sub reallign {
	(my $WORKING_DIR,
	 my $SAMPLE_NAME,
	 my $RefGenome,
	 my $ALIGNMENT_FP)
	= @_;

	# find Picard
	my 	$PicardJarPath = `which picard`;
		$PicardJarPath =~ s/\n//g;
		$PicardJarPath =~ s/bin\/.*$/bin\//g;
		$PicardJarPath .= "../share/picard-2.20.8-0/picard.jar";

	# create folder
	unless(-d "$WORKING_DIR/data/"){system("mkdir $WORKING_DIR/data/");}

	if ( ! -f $ALIGNMENT_FP ){
		print "Error the bam wasn't find";
		die;
	}
	# bam to fastq
	system("java -Xmx8G -jar $PicardJarPath SamToFastq I=$ALIGNMENT_FP FASTQ=$WORKING_DIR/data/read_1.fq SECOND_END_FASTQ=$WORKING_DIR/data/read_2.fq");
	# align on new reference
	system("bwa mem $RefGenome  $WORKING_DIR/data/read_1.fq $WORKING_DIR/data/read_2.fq | samtools sort  -o $WORKING_DIR/data/align.bam");

	# to avoid an error
	system("java  -Xmx8G -jar $PicardJarPath ReplaceSamHeader I= $ALIGNMENT_FP HEADER=$WORKING_DIR/REF/hg38.chr.dict O=$WORKING_DIR/data/$SAMPLE_NAME.align.bam");
	# indexing
	system("samtools index $WORKING_DIR/data/$SAMPLE_NAME.align.bam");
}

sub GenerateMaskedAlignmentAndVcf {

	# recuperate the variable
	(my $WORKING_DIR,
	 my $SAMPLE_NAME,
	 my $ALIGNMENT_FP,
	 my $NR_OF_THREADS)
	= @_;
	
	# find Picard
	my 	$PicardJarPath = `which picard`;
		$PicardJarPath =~ s/\n//g;
		$PicardJarPath =~ s/bin\/.*$/bin\//g;
		$PicardJarPath .= "../share/picard-2.20.8-0/picard.jar";

	# create the variable
	my 	$ExtractionBedFP		= "";
	my 	$VarCallBedFP			= "";
	my 	$CovExonsFP				= "";
	my 	$RefGenome				= "";
	my 	$MaskedRefGenome		= "";
	my 	$GenomeFileFilePath 	= "";
	
	# test if working directory, refseq and BEDs exists
	if (! -d $WORKING_DIR){
		print "$WORKING_DIR should already exist.\n";
		print "Please enter an existing working directory.\n";
		die("$!\n");
	}

	# Recuperate the path of all files

	$ExtractionBedFP		= "$WORKING_DIR/BED/all.forextraction.hg38.bed";
	$VarCallBedFP			= "$WORKING_DIR/BED/all.forvarcall.hg38.bed";
	$CovExonsFP				= "$WORKING_DIR/BED/all.homologousexons.hg38.bed";
	$RefGenome				= "$WORKING_DIR/REF/hg38.chr.fa";
	$MaskedRefGenome		= "$WORKING_DIR/REF/hg38.masked.chr.fa";
	$GenomeFileFilePath		= "$WORKING_DIR/REF/hg38.chr.txt";


	# if one files is missing
	if ( ! -f $ExtractionBedFP  	||
		 ! -f $VarCallBedFP  		||
		 ! -f $CovExonsFP  			||
		 ! -f $RefGenome 			||
		 ! -f $MaskedRefGenome		||
		 ! -f $GenomeFileFilePath){
		
		print "BEDs and RefSeqs should already be present\n";
		print "Please first run PrepareBED and MaskReferenceGenome\n";
		die("$!\n");	
	}

	# Create the folder
	unless(-d "$WORKING_DIR/RAW/"){system("mkdir $WORKING_DIR/RAW/");}
	unless(-d "$WORKING_DIR/RAW/"){system("mkdir $WORKING_DIR/RAW/$SAMPLE_NAME");}

	### First part create vcf files

	# alignement to bam (taking the reads and genome) positin? name reads
	system	("samtools view -b -L $ExtractionBedFP -o $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.ori.bam $ALIGNMENT_FP");
	# classified with the coordinated
	system	("java -Xmx8G -jar $PicardJarPath SortSam I=$WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.ori.bam O=$WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.ori.sorted.bam SO=coordinate");
	# tagg the duplicated from the files
	system	("java -Xmx8G -jar $PicardJarPath MarkDuplicates I=$WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.ori.sorted.bam O=$WORKING_DIR/$SAMPLE_NAME/$SAMPLE_NAME.ori.sorted.remdup.bam M=$WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.ori.sorted.remdup.met REMOVE_DUPLICATES=true");
	# remove the sorted bam file
	system	("rm -f $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.ori.sorted.bam");
	# bam to bai
	system	("samtools index $WORKING_DIR/$SAMPLE_NAME/$SAMPLE_NAME.ori.sorted.remdup.bam");
	# eliminate the reads to have uniq reads
	system	("samtools view -bh -q 1 $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.ori.sorted.remdup.bam > $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.ori.sorted.remdup.uniq.bam");
	system("echo second");
	# compare the exons homologous with the reads and counts -> put in a zip
	system	("bedtools coverage -a $CovExonsFP -b  $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.ori.sorted.remdup.bam -counts | bgzip -c >  $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.ori.sorted.remdup.cov.gz");
	# compare the exons homologous  with the uniq reads and counts -> put in a zip
	system	("bedtools coverage -a $CovExonsFP -b  $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.ori.sorted.remdup.uniq.bam -counts | bgzip -c >  $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.ori.sorted.remdup.uniq.cov.gz");
	# Infer the SNV and indel from reads (here)
	system("echo third");
	system	("lofreq call --call-indels --force-overwrite --no-default-filter --use-orphan -a 1 -b 20 -B -N -f $RefGenome -l $VarCallBedFP -o  $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.ori.sorted.remdup.lofreq.vcf  $WORKING_DIR/RAW/$SAMPLE_NAME.ori.sorted.remdup.bam");

	system("echo third");
	# zip the inference 
	system	("bgzip -c  $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.ori.sorted.remdup.lofreq.vcf >  $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.ori.sorted.remdup.lofreq.vcf.gz");
	# Index and retirve overlapping specified regions
	system	("tabix -p vcf  $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.ori.sorted.remdup.lofreq.vcf.gz");
	# remove the inference file
	system	("rm -f  $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.ori.sorted.remdup.lofreq.vcf");
	# deteckt the oitential variant sites per samples
	system("java -jar $PicardJarPath AddOrReplaceReadGroups I=$WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.ori.sorted.remdup.bam O=$WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.ori.asorted.remdup.bam RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=hello");
	system("samtools index $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.ori.asorted.remdup.bam");
	system	("gatk --java-options \-Xmx4g\ HaplotypeCaller -R $RefGenome -I  $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.ori.asorted.remdup.bam -O  $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.ori.sorted.remdup.gatk.vcf -L $VarCallBedFP --annotation MappingQualityRankSumTest --annotation MappingQualityZero --annotation QualByDepth --annotation ReadPosRankSumTest --annotation RMSMappingQuality --annotation BaseQualityRankSumTest --annotation FisherStrand --annotation MappingQuality --annotation DepthPerAlleleBySample");
	system("echo yes");

	### Clean the files 
	# zip the results
	system	("bgzip -c  $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME\.ori.sorted.remdup.gatk.vcf >  $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.ori.sorted.remdup.gatk.vcf.gz");
	# Index and retirve overlapping specified regions
	system	("tabix -p vcf  $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME\.ori.sorted.remdup.gatk.vcf.gz");
	# removes files
	system	("rm -f  $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.ori.sorted.remdup.gatk.vcf");
	system	("rm -f  $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.ori.sorted.remdup.gatk.vcf.idx");

	###Second part
	system("echo fourth");
	# compare the VarCallBedFP with the uniq reads and counts -> put in a zip
	system	("bedtools coverage -a $VarCallBedFP -b  $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.ori.sorted.remdup.bam -d -sorted  -g $GenomeFileFilePath >  $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.ori.sorted.remdup.cov");
	system	("bgzip -c  $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.ori.sorted.remdup.cov >  $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.ori.sorted.remdup.cov.gz");
	system	("tabix -p bed  $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.ori.sorted.remdup.cov.gz");
	system	("rm -f  $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.ori.sorted.remdup.cov");
	
	# goups read together by names and the way of reads
	system	("samtools collate -uO  $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.ori.sorted.remdup.bam  $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME | samtools fastq -1  $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.1.fastq -2  $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.2.fastq -0  $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.0.fastq -t -");
	system	("gzip  $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.1.fastq");
	system	("gzip  $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.2.fastq");
	# Re-pairs reads that became disordered or had some mates eliminated.
	system	("repair.sh -in1=$WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.1.fastq.gz in2=$WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.2.fastq.gz out1=$WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.rep.1.fastq.gz out2=$WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.rep.2.fastq.gz outsingle=$WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.rep.0.fastq.gz usejni=t");
	system	("rm -f  $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.1.fastq.gz");
	system	("rm -f  $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.2.fastq.gz");
	system	("rm -f  $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.0.fastq");

	system("echo fifth");
	# Local aligment 	
	system	("bwa mem -t $NR_OF_THREADS -M $MaskedRefGenome  $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.rep.1.fastq.gz  $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.rep.2.fastq.gz >  $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.masked.sam");
	# S for compatibilty not needed b -> bam format
	system	("samtools view -Sb  $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.masked.sam >  $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.masked.bam");
	# Class with coordinate
	system	("java -Xmx8G -jar $PicardJarPath SortSam I= $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.masked.bam O= $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.masked.sorted.bam SO=coordinate");
	
	system	("rm -f $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.rep.0.fastq.gz");
	system	("rm -f $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.rep.1.fastq.gz");
	system	("rm -f $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.rep.2.fastq.gz");
	system	("rm -f $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.masked.sam");
	system	("rm -f $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.masked.bam");

	# Erase duplicate from file 
	system	("java -Xmx8G -jar $PicardJarPath MarkDuplicates I=$WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.masked.sorted.bam O=$WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.masked.sorted.remdup.bam M=$WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.masked.sorted.remdup.met REMOVE_DUPLICATES=true");
	system	("rm -f $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.masked.sorted.bam");
	# add indel quality 
	system	("lofreq  indelqual --dindel -f $MaskedRefGenome -o $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.masked.sorted.remdup.bqsr.bam $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.masked.sorted.remdup.bam");
	
	system	("rm -f $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.masked.sorted.remdup.bam");

	# index 
	system	("samtools index $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.masked.sorted.remdup.bqsr.bam");
	# Assigns all the reads in a file to a single new read-group
	system	("java -Xmx8G -jar $PicardJarPath AddOrReplaceReadGroups I=$WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.masked.sorted.remdup.bqsr.bam O=$WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.masked.sorted.remdup.bqsr.rg.bam RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$SAMPLE_NAME");			
	system("echo sixth");
	system	("rm -f $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.masked.sorted.remdup.bqsr.bam");
	system	("rm -f $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.masked.sorted.remdup.bqsr.bam.bai");
	# index 
	system	("samtools index $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.masked.sorted.remdup.bqsr.rg.bam");
	### third part
	# Infer the SNV and indel from reads (here)
	system	("lofreq call --call-indels --force-overwrite --no-default-filter --use-orphan -a 1 -b 20 -B -N -f $MaskedRefGenome -l $VarCallBedFP -o $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.masked.sorted.remdup.bqsr.rg.lofreq.vcf $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.masked.sorted.remdup.bqsr.rg.bam");
	system	("bgzip -c $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.masked.sorted.remdup.bqsr.rg.lofreq.vcf > $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.masked.sorted.remdup.bqsr.rg.lofreq.vcf.gz");
	system	("tabix -p vcf $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.masked.sorted.remdup.bqsr.rg.lofreq.vcf.gz");
	system	("rm -f $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.masked.sorted.remdup.bqsr.rg.lofreq.vcf");
	#  In other words, whenever the program encounters a region showing signs of variation, it discards the existing mapping information and completely reassembles the reads in that region. 
	system	("gatk --java-options \-Xmx4g\ HaplotypeCaller -R $MaskedRefGenome -I $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.masked.sorted.remdup.bqsr.rg.bam -O $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.masked.sorted.remdup.bqsr.rg.gatk.vcf -L $VarCallBedFP --annotation MappingQualityRankSumTest --annotation MappingQualityZero --annotation QualByDepth --annotation ReadPosRankSumTest --annotation RMSMappingQuality --annotation BaseQualityRankSumTest --annotation FisherStrand --annotation MappingQuality --annotation DepthPerAlleleBySample");
	
	system	("bgzip -c $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.masked.sorted.remdup.bqsr.rg.gatk.vcf > $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.masked.sorted.remdup.bqsr.rg.gatk.vcf.gz");
	system	("tabix -p vcf $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.masked.sorted.remdup.bqsr.rg.gatk.vcf.gz");
	system	("rm -f $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.masked.sorted.remdup.bqsr.rg.gatk.vcf");
	system	("rm -f $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.masked.sorted.remdup.bqsr.rg.gatk.vcf.idx");	
	system	("bedtools coverage -a $VarCallBedFP -b $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.masked.sorted.remdup.bqsr.rg.bam -d -sorted  -g $GenomeFileFilePath > $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.masked.sorted.remdup.bqsr.rg.cov");
	
	system	("bgzip -c $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.masked.sorted.remdup.bqsr.rg.cov > $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.masked.sorted.remdup.bqsr.rg.cov.gz");
	system	("tabix -p bed $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.masked.sorted.remdup.bqsr.rg.cov.gz");
	system	("rm -f $WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.masked.sorted.remdup.bqsr.rg.cov");
}

sub FilterRawVariants {
	
	(my $WORKING_DIR,
	 my $SAMPLE_NAME)
	= @_;
	
	my $SDsBedFP						= "";
	my $AllToMapOnToOtherFP				= "";
	my $MappedSDsFP						= "";
	my $PositionToRegionFP				= "";
	my $HPsBedFP						= "";
	my $TroubleSitesFP					= "";
	my $CohortAF_FP						= "";
	
	my $RegionToStrandsFP				= "$WORKING_DIR/BED/RegionID_ToStrand.txt";		
	my $VariantsOutFP					= "$WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.RawSNVs.txt";
	my $VariantsFilter1OutFP			= "$WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.Filter1SNVs.txt";
	my $VariantsFilter2OutFP			= "$WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.Filter2SNVs.txt";
	my $VariantsFilter1BedFP			= "$WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.Filter1SNVs.bed";
	my $SortedVariantsFilter1BedFP		= "$WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.Filter1SNVs.sorted.bed";
	my $VariantsFilter1HpAnnoBedFP		= "$WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.Filter1SNVs.HpAnno.bed";
	my $FileHandle						= new IO::Zlib;
	my %SDsBED 							= ();
	my %SDsVAR 							= ();
	my @VCFs							= ( "$SAMPLE_NAME.masked.sorted.remdup.bqsr.rg.lofreq.vcf.gz",
											"$SAMPLE_NAME.ori.sorted.remdup.gatk.vcf.gz");
	my $CovFilePath						= "$WORKING_DIR/RAW/$SAMPLE_NAME/$SAMPLE_NAME.masked.sorted.remdup.bqsr.rg.cov.gz";
	my %LoFreq_Variants					= ();
	my %GATK_Variants					= ();
	my %PosOfInterest					= ();
	my %CoordinateMappings				= ();
	my %Coverage						= ();
	my %MappedSDs						= ();
	my %PosToRegionName					= ();
	my %Strands							= ();
	my %VarsToExclude					= ();
	my %F1Variants						= ();
	my %SDtroublePos					= ();
	my %PopFreq							= ();

	$SDsBedFP						= "$WORKING_DIR/BED/AllSDs.noalt.chr.bed";
	$AllToMapOnToOtherFP			= "$WORKING_DIR/BED/AllToMapOnToOther.chr.txt.gz";	
	$MappedSDsFP					= "$WORKING_DIR/BED/MappedSDs.chr.txt";
	$PositionToRegionFP				= "$WORKING_DIR/BED/PosToRegionID.chr.txt.gz";
		
	$HPsBedFP						= "$WORKING_DIR/BED/HPs.noalt.chr.bed";
	$TroubleSitesFP					= "$WORKING_DIR/BED/SitesToExclude.noalt.chr.bed";
	$CohortAF_FP					= "$WORKING_DIR/BED/CohortAlleleFreq.chr.txt";


	# Read in SDs #
	
	print "Read in SDs\n";
	# Prepare
	open SD, "$SDsBedFP" or die ("Can't open $SDsBedFP\n");
	while(<SD>){
		# Erase the newline
		$_ =~ s/\n//g;
		# Erase the cariage
		$_ =~ s/\r//g;
			
		(my $Chrom, undef, undef, my $SD) = split (m/\t/, $_);
		# chr1~879153~879153~C/T
		(undef, my $VarStart, undef, undef) = split(m/~/, $SD);

		$SDsBED{$Chrom}{$VarStart}{$SD} = undef;
		$SDsVAR{$SD} = undef;		
		# chrom posstart SD = undef
		# SD = undef
	}
	close SD;
	
	# Read In coorindate mappings #
	
	print "Read In coorindate mappings\n";
	
	$FileHandle	= new IO::Zlib;
	
	if ($FileHandle->open($AllToMapOnToOtherFP, "rb")){
		# chr10~102034807	chrX~63946242
		while (<$FileHandle>){
			# Erase the newline
			$_ =~ s/\n//g;
			# chr10~102034807	chrX~63946242	
			(my $ToMapOnCoord, my $OtherCoord) = split(m/\t/, $_);
		
			$CoordinateMappings{$ToMapOnCoord}{$OtherCoord} = undef;	
			# CoordinateMappings start end = undef
		}
	}
	
	# Read in mapped SDs
	
	print "Read in mapped SDs\n";
	
	open M, "$MappedSDsFP" or die ("Can't open $MappedSDsFP\n");
	while (<M>){
	
		$_ =~ s/\n//g;
		#chr10~102034829~102034829~C/-	chrX~63946264~63946263~-/C
		(my $ToMapOnSD, my $OtherSD) = split (m/\t/, $_);
		
		$MappedSDs{$ToMapOnSD}{$OtherSD} = undef;
		#chr10~102034829~102034829~C/-	chrX~63946264~63946263~-/C = undef
	}
	close M;
	
	# Read in position to region ID
	
	print "Read in position to region ID\n";
	
	$FileHandle	= new IO::Zlib;

	if ($FileHandle->open($PositionToRegionFP, "rb")){
			
		while (<$FileHandle>){
								
			$_ =~ s/\n//g;
			#chr1~877479	SAMD11
			(my $Coord, my $RegionName) = split(m/\t/, $_);
			#chr1            877479	
			(my $Chrom, my $Position) = split(m/~/, $Coord);
		
			$PosToRegionName{$Chrom}{$Position} = $RegionName;	
			# chr1 877479 = SAMD11
		}
	}
	
	# Read in strands
	
	print "Read in strands\n";
	
	open S, "$RegionToStrandsFP" or die ("Can't open $RegionToStrandsFP\n");
	while(<S>){
		
		$_ =~ s/\n//g;
		# PRDX2P2	+
		(my $RegionID, my $Strand) = split(m/\t/, $_);
		
		$Strands{$RegionID} = $Strand;
		# PRDX2P2 =	+
	}
	close S;
	
	# Tabulate variants #
	
	print "Tabulate variants\n";
	# Open raw SNV
	open O, ">$VariantsOutFP" or die ("Can't open $VariantsOutFP\n");
	foreach my $VCF (@VCFs){
		
		my 	$Caller 						= "";
		my 	$Masking						= "";
		my 	$Analysis						= "";
		my 	$VARIANT_CALLER_SELECT			= "";
			$FileHandle						= new IO::Zlib;
		
		my $VCF_FilePath 					= "$WORKING_DIR/RAW/$SAMPLE_NAME/$VCF";
		
		if 		($VCF =~ m/\.gatk\./)	{$Caller = "gatk";}
		else 							{$Caller = "lofreq";}
		
		if 		($VCF =~ m/\.masked\./)	{$Analysis = "MASKED";}
		elsif 	($VCF =~ m/\.ori\./)	{$Analysis = "ORI";}
		
		if 		($Caller eq "lofreq")	{$VARIANT_CALLER_SELECT = "LoFreq";}
		elsif 	($Caller eq "gatk")		{$VARIANT_CALLER_SELECT = "GATK";}
		
		if ($FileHandle->open($VCF_FilePath, "rb")){
			
			while (<$FileHandle>){
			
				next if ($_ =~ m/^#/);
									
				my 	$Line 	= $_;
					$Line 	=~ s/\n//g;
				# here a function to create the table of variabt chr1~879153~879153~C/T 
				(my $VariantsRef, my $AFsRef, my $Depth, my $Qual) = DetermineVariantInVCF_Line ($Line, $VARIANT_CALLER_SELECT);
				
				for (my $I = 0; $I < scalar @$VariantsRef; $I++){
			
					my 	$Variant = $VariantsRef->[$I];
													
					(my $Chrom, my $ChromStart, my $ChromEnd, undef) = split(m/~/, $Variant);
					
					my 	$AF				= $AFsRef->[$I];
						$AF 			=~ s/AF\=//g;
					my 	$Min			= min($ChromStart, $ChromEnd);
					my 	$Max			= max($ChromStart, $ChromEnd);
					my 	$BedStart 		= $Min-1;
					my 	$VariantClass	= "nonSD";
					# look if variant exist in SD
					if 		(exists $SDsVAR{$Variant}){$VariantClass = "SD";}
					
					next if ($AF eq "/");
													
					$AF 								 = sprintf("%.3f", $AF);
					
					print O "$Chrom\t$BedStart\t$Max\t$Variant\t$SAMPLE_NAME\t$Depth\t$AF\t$Qual\t$Caller\t$Analysis\t$VariantClass\n";
					
					if  ($Caller eq "lofreq"){
						
						$LoFreq_Variants{$Chrom}{$ChromStart}{$Variant}{$Analysis} = $AF;
					}
					
					if  ($Caller eq "gatk"){
						
						$GATK_Variants{$Chrom}{$ChromStart}{$Variant}{$Analysis} = $AF;
					}	
				}
			}
		}
		$FileHandle->close;
	}
	close O;
	
	# Positions of interest obtain in tabulate variants#
	
	print "Positions of interest\n";
	# from file lofreq
	foreach my $Chrom (keys %LoFreq_Variants){
		foreach my $Start (keys %{$LoFreq_Variants{$Chrom}}){
			foreach my $Variant (keys %{$LoFreq_Variants{$Chrom}{$Start}}){
				
				$PosOfInterest{$Chrom}{$Start}{$Variant} = undef;
			}
		}
	}
	# from file gatk
	foreach my $Chrom (keys %GATK_Variants){
		foreach my $Start (keys %{$GATK_Variants{$Chrom}}){
			foreach my $Variant (keys %{$GATK_Variants{$Chrom}{$Start}}){
				
				$PosOfInterest{$Chrom}{$Start}{$Variant} = undef;
			}
		}
	}
	# from file SG
	foreach my $Chrom (keys %SDsBED){
		foreach my $Start (keys %{$SDsBED{$Chrom}}){
			foreach my $Variant (keys %{$SDsBED{$Chrom}{$Start}}){
				
				$PosOfInterest{$Chrom}{$Start}{$Variant} = undef;
			}
		}
	}
	
	# Correct with coverage above threshold #
	
	print "Correct with coverage above threshold\n";
	
	$FileHandle	= new IO::Zlib;
	
	if ($FileHandle->open($CovFilePath, "rb")){
			
		while (<$FileHandle>){
								
			my 	$Line 	= $_;
				$Line 	=~ s/\n//g;
			
			(my $Chrom, my $Start, my $End, my $PosInInt, my $Cov) = split (m/\t/, $Line);
			
			my $CovPos = $Start + $PosInInt;
			# if not in interresting regions
			next if (!exists $CoordinateMappings{"$Chrom\~$CovPos"});
			# quality enough
			if ($Cov >= 60){

				$Coverage{$Chrom}{$CovPos} = $Cov; 
							
				if (	exists $LoFreq_Variants{$Chrom}						&&
						exists $LoFreq_Variants{$Chrom}{$CovPos}){
							
					foreach my $Variant (keys %{$LoFreq_Variants{$Chrom}{$CovPos}}){
						# initialize masked if doesn't exist in lofreq
						if (!exists $LoFreq_Variants{$Chrom}{$CovPos}{$Variant}{"masked"}){
						
							$LoFreq_Variants{$Chrom}{$CovPos}{$Variant}{"masked"} = 0;
						}
					}
					
					foreach my $Variant (keys %{$PosOfInterest{$Chrom}{$CovPos}}){
						# initialize masked if doesn't exist in position of interesst
						if (!exists $LoFreq_Variants{$Chrom}{$CovPos}{$Variant} &&
							!exists $LoFreq_Variants{$Chrom}{$CovPos}{$Variant}{"masked"}){
							
							$LoFreq_Variants{$Chrom}{$CovPos}{$Variant}{"masked"} = 0;
						}
					}						
				}
				elsif 	((!exists $LoFreq_Variants{$Chrom}						||
						  !exists $LoFreq_Variants{$Chrom}{$CovPos}) 					&&
						  exists $PosOfInterest{$Chrom}									&&
						  exists $PosOfInterest{$Chrom}{$CovPos}){
					# initialized masked if doesn't exist in lofreq but exist in pos interrested
					foreach my $Variant (keys %{$PosOfInterest{$Chrom}{$CovPos}}){
													
						$LoFreq_Variants{$Chrom}{$CovPos}{$Variant}{"masked"} = 0;
					}
				}
			}
		}
	}
	$FileHandle->close;	
	
	# Correct with coverage below threshold #
	
	print "Correct with coverage below threshold\n";
	
	foreach my $Chrom (sort keys %PosOfInterest){
		foreach my $Start (sort keys %{$PosOfInterest{$Chrom}}){
		
			next if (!exists $CoordinateMappings{"$Chrom\~$Start"});
			# exist if in lofreq
			if 	(	exists $LoFreq_Variants{$Chrom} &&
					exists $LoFreq_Variants{$Chrom}{$Start}){
				# delete the lofreq variant if the cov quality under 60 
				if (!exists $Coverage{$Chrom}{$Start}){
					
					delete($LoFreq_Variants{$Chrom}{$Start});
					
					foreach my $HomoloCoord (keys %{$CoordinateMappings{"$Chrom\~$Start"}}){
						
						(my $HomChrom, my $HomPos) = split(m/~/, $HomoloCoord);
						# delete the lofreq variant if present in coordinate mapping
						if ($LoFreq_Variants{$HomChrom}{$HomPos}){
								
							delete($LoFreq_Variants{$HomChrom}{$HomPos});
						}
					}
				}
			}
		}
	}
		
	open F1, ">$VariantsFilter1OutFP" or die ("can't open $VariantsFilter1OutFP\n");
	open B, ">$VariantsFilter1BedFP" or die ("can't open $VariantsFilter1BedFP\n");
	print F1 "Chromosome\tStart\tVAP\tVariantCall\tVAF_Masked\tMaskedCov\n";

	foreach my $Chrom (sort keys %LoFreq_Variants){
		foreach my $Pos (sort keys %{$LoFreq_Variants{$Chrom}}){
			
			#print $Pos . "\n";
			# pass if quality high
			next if (!exists $Coverage{$Chrom}{$Pos});
			
			foreach my $Variant (sort keys %{$LoFreq_Variants{$Chrom}{$Pos}}){
				
				#print $Variant . "\n";
				#
				(my $ToMapOnChrom, my $ToMapOnStart, my $ToMapOnEnd, my $Change) = split (m/~/, $Variant);
				# mutation
				(my $RefNucl, my $AltNucl) 		= split(m/\//, $Change);
				# if have pos start and pos end
				if ($CoordinateMappings{"$Chrom\~$ToMapOnStart"} && $CoordinateMappings{"$Chrom\~$ToMapOnEnd"}){	
					# exist in SD; in lofreq ;qulity superior at 0.15 but doesnt exist in GATK initlilize present and edge
					if (! exists $MappedSDs{$Variant}															&&
						  exists $LoFreq_Variants{$Chrom}{$Pos}{$Variant}{"MASKED"} 					&&
								 $LoFreq_Variants{$Chrom}{$Pos}{$Variant}{"MASKED"} >=0.15				&&
						(! exists $GATK_Variants{$Chrom} 												||
						 ! exists $GATK_Variants{$Chrom}{$Pos}											||
						 ! exists $GATK_Variants{$Chrom}{$Pos}											||
						 ! exists $GATK_Variants{$Chrom}{$Pos}{$Variant}								||
						 ! exists $GATK_Variants{$Chrom}{$Pos}{$Variant}{"ORI"})){
						
						my $Present				= 0;
						my $Edge				= 0;
						
						### CHECK NEIGHBOURHOOD FOR INDELS ###
						
						if ($RefNucl eq "-" || $AltNucl eq "-"){
										
							for (my $I = $Pos-2; $I <= $Pos+2; $I++){
								# exist variant near mut check 
								if (exists $GATK_Variants{$Chrom} 						&&
									exists $GATK_Variants{$Chrom}{$I}){
									
									foreach my $TryVariant (keys %{$GATK_Variants{$Chrom}{$I}}){
										
										(undef, undef, undef, my $TryChange) 	= split(m/\~/, $TryVariant);
										(my $TryRefNucl, my $TryAltNucl) 		= split(m/\//, $TryChange);
										if (exists $GATK_Variants{$Chrom}{$I}{$TryVariant}{"ORI"} &&
											($RefNucl eq "-" && $TryRefNucl eq "-") ||
											($AltNucl eq "-" && $TryAltNucl eq "-")){
											
											$Present = 1;
										}
									}
								}
							}
						}
							
						## CHECK IF ALL OTHERVARIANTS ARE ABSENT ###
								
						my $OtherVariantsRef = ConstructOtherVariant($Variant, \%CoordinateMappings, \%PosToRegionName, \%Strands, \%SDsBED);
						# flag 
						if ($Present == 0){
							
							foreach my $OtherVariant (keys %$OtherVariantsRef){
							
								(my $OtherChrom, my $OtherStart, my $OtherEnd, my $OtherChange) = split(m/\~/, $OtherVariant);
								# exist gatk variant
								if (exists $GATK_Variants{$OtherChrom} 											&&
									exists $GATK_Variants{$OtherChrom}{$OtherStart}								&&
									exists $GATK_Variants{$OtherChrom}{$OtherStart}{$OtherVariant}				&&
									exists $GATK_Variants{$OtherChrom}{$OtherStart}{$OtherVariant}{"ORI"}){
									
									$Present = 1;
								}
								# doesn t exist postoregion then flag
								if (!exists $PosToRegionName{$OtherChrom} ||
									!exists $PosToRegionName{$OtherChrom}{$OtherStart}){
									
									$Edge = 1;
								}
								
								if ($Present+$Edge == 0 && ($RefNucl eq "-" || $AltNucl eq "-")){
									
									for (my $I = $OtherStart-2; $I <= $OtherStart+2; $I++){
										
										if (exists $GATK_Variants{$OtherChrom} 					&&
											exists $GATK_Variants{$OtherChrom}{$I}					&&
											exists $GATK_Variants{$OtherChrom}{$I}){
											
											foreach my $TryVariant (keys %{$GATK_Variants{$OtherChrom}{$I}}){
												
												(undef, undef, undef, my $TryChange) 	= split(m/\~/, $TryVariant);
												(my $TryRefNucl, my $TryAltNucl) 		= split(m/\//, $TryChange);
												if (exists $GATK_Variants{$OtherChrom}{$I}{$TryVariant}{"ORI"} &&
													($RefNucl eq "-" && $TryRefNucl eq "-") ||
													($AltNucl eq "-" && $TryAltNucl eq "-")){
													
													$Present = 1;
												}
											}
										}
									}
								}
							}
						}
							
						if ($Present+$Edge == 0){
							
							print F1 $Chrom . "\t" . $Pos . "\t" . $Variant . "\t" . $Variant . "\t" . $LoFreq_Variants{$Chrom}{$Pos}{$Variant}{"MASKED"} . "\t" . $Coverage{$Chrom}{$Pos} . "\n";
							
							(undef, my $StartPos, my $EndPos, undef) = split(m/\~/, $Variant);
		
							my $BedStart = min($StartPos, $EndPos) - 1;
							my $BedEnd = min($StartPos, $EndPos);
		
							print B $Chrom . "\t" . $BedStart . "\t" . $BedEnd . "\t" . $Variant . "\n"; 
							
							$F1Variants{$Variant} = "$Chrom\t$Pos\t$Variant\t$Variant\t$LoFreq_Variants{$Chrom}{$Pos}{$Variant}{\"MASKED\"}\t$Coverage{$Chrom}{$Pos}";
							
							foreach my $OtherVariant (keys %$OtherVariantsRef){
						
								(my $OtherChrom, my $OtherStart, my $OtherEnd, my $OtherChange) = split(m/\~/, $OtherVariant);
								
								my $OtherRegionName = $PosToRegionName{$OtherChrom}{$OtherStart};
								
								print F1 $OtherChrom . "\t" . $OtherStart . "\t" . $OtherVariant . "\t" . $Variant . "\t" . $LoFreq_Variants{$Chrom}{$Pos}{$Variant}{"MASKED"}  . "\t" . $Coverage{$Chrom}{$Pos} . "\n";	
								
								my $OtherBedStart = min($OtherStart, $OtherEnd) - 1;
								my $OtherBedEnd = min($OtherStart, $OtherEnd);
		
								print B $OtherChrom . "\t" . $OtherBedStart . "\t" . $OtherBedEnd . "\t" . $OtherVariant . "\n"; 
								
								$F1Variants{$OtherVariant} = "$OtherChrom\t$OtherStart\t$OtherVariant\t$Variant\t$LoFreq_Variants{$Chrom}{$Pos}{$Variant}{\"MASKED\"}\t$Coverage{$Chrom}{$Pos}";
							}	
						}
					}					
				}
			}
		}	
	}
	
	close B;
	close F1;
	
	
	# Filter Out HP regions #
	
	system("sort-bed $VariantsFilter1BedFP > $SortedVariantsFilter1BedFP");
	system("bedmap --echo --echo-map --delim '\t' $SortedVariantsFilter1BedFP $HPsBedFP > $VariantsFilter1HpAnnoBedFP");
	
	open HP, "$VariantsFilter1HpAnnoBedFP" or die ("Can't open $VariantsFilter1HpAnnoBedFP\n");
	while(<HP>){
		
		my 	$Line 				= $_;
			$Line 				=~ s/\n//g;
		my 	@LineValues			= split(m/\t/, $Line);
		my 	$Variation			= $LineValues[3];
		
		if (scalar @LineValues == 8){
						
			$VarsToExclude{$Variation} = undef;
		}
	}
	close HP;
	
	# Trouble regions #
	
	open S, "$TroubleSitesFP" or die ("Can't open $TroubleSitesFP\n");
	while(<S>){
		
		$_ =~ s/\n//g;
		
		(my $Chrom, undef, my $Pos) = split(m/\t/, $_);
		
		$SDtroublePos{$Chrom}{$Pos} = undef;	
	}
	close S;
	
	# Population frequency #
	
	open P, "$CohortAF_FP" or die ("Can't open $CohortAF_FP");
	while(<P>){
		
		$_ =~ s/\n//g;
		
		(my $Variation, my $Frequency) = split(m/\t/, $_);
		
		$PopFreq{$Variation} = $Frequency;
	}
	close P;
	
	# Write out F2 variants #
	
	open F2, ">$VariantsFilter2OutFP" or die ("can't open $VariantsFilter2OutFP\n");
	print F2 "Chromosome\tStart\tVAP\tVariantCall\tVAF_Masked\tMaskedCov\tPopFreq\n";
	
	foreach my $Variant (sort keys %F1Variants){
		
		if (! exists $VarsToExclude{$Variant}){
			
			(my $Chrom, my $Start) = split(m/\~/, $Variant);
			
			if (! exists $SDtroublePos{$Chrom}{$Start}){
				
				if (exists $PopFreq{$Variant}){
					print F2 $F1Variants{$Variant} . "\t" . $PopFreq{$Variant} . "\n";
				}
				else {
					print F2 $F1Variants{$Variant} . "\t" . "0" . "\n";
				}
				
				
			}
		}
	}
	
	close F2;
}
	
	
	




# determined the type of mutations
sub DetermineVariantInVCF_Line {
	# line of the file; which type of file
	(my $Line,
	 my $Caller)
	= @_;
	
	my @LineValues 			= split ("\t", $Line);
	
	my $Chrom				= $LineValues[0];
	my $Pos					= $LineValues[1];
	my $ReferenceBase		= $LineValues[3];
	my $VariantBases		= $LineValues[4];
	my $Qual				= $LineValues[5];
	my $Info 				= $LineValues[7];
	my $Depth				= "/";
	my $Variant				= "/";
	
	my @AFs					= ();
	my @Variants			= ();
	my @VariantBases		= split (m/\,/, $VariantBases);
	
	foreach my $VariantBase (@VariantBases){
	
		### Determine Variant ###

		if (length($ReferenceBase) == length ($VariantBase)){ 		# SUBSTITUTIONS
			
			$Variant = $Chrom . "~" . $Pos . "~" . $Pos . "~" .  "$ReferenceBase\/$VariantBase";
		}
		elsif (length($ReferenceBase) < length ($VariantBase)){ 	# INSERTION
			
			($ReferenceBase, $VariantBase, $Pos) = FormatGATKCoordinatesOfVariant($ReferenceBase, $VariantBase, $Pos);
			
			my $StartPos 	= $Pos;
			my $EndPos  	= $Pos - 1;
			
			$Variant = $Chrom . "~" . $StartPos . "~" . $EndPos . "~" .  "$ReferenceBase\/$VariantBase";
		}
		else { 														# DELETION
		
			($ReferenceBase, $VariantBase, $Pos) = FormatGATKCoordinatesOfVariant($ReferenceBase, $VariantBase, $Pos);
			
			my $StartPos = $Pos;
			my $EndPos   = $Pos + length($ReferenceBase) - 1;
			
			$Variant = $Chrom . "~" . $StartPos . "~" . $EndPos . "~" .  "$ReferenceBase\/$VariantBase";
		}
				
		push (@Variants, $Variant);
	}
	
	##########################
	
	if 		($Caller =~	m/GATK/){
		
		my 	$Format 	= $LineValues[8];
		my 	$Sample 	= $LineValues[9];
		
		# print $Caller . "\t" . $Format . "\t" . $Sample . "\n";
		
		if (defined $Format && $Format eq "GT:AD:DP:GQ:PL"){
		
			$Sample =~ m/^.*:(.*):(.*):.*:.*$/;
			
			if (!$1){
			
				print $Caller . "\t" . $Line . "\n";
			}
			
			# print $1 . \t . $2 . \n;
			
			my 	@ADs 	= split (m/\,/, $1);
				$Depth	= $2;
			
			my $RefAD  	= shift (@ADs);
			
			foreach my $AltAD (@ADs){
			
				if ($AltAD+$RefAD==0){
					push(@AFs, "/");
				}
				else {
					push(@AFs, ($AltAD/($AltAD+$RefAD)));
				}
			}
		}
		elsif ($LineValues[7] =~ m/\;AF\=.+?\;.*/g){
		
			$LineValues[7] =~ m/\;AF=(.+?)\;/;
			
			@AFs = split (m/\,/, $1);
		}
		else {
		
			if ($Line =~ m/\tPASS\t/){
				print $Line . "\n";
				
				die ("Format does not match: GT:AD:DP:GQ:PL\n");
			}
		}							
	}
	elsif 	($Caller eq "LoFreq"){
		
		($Depth, my $AFs, undef, undef) = split (/\;/, $Info);
		$Depth =~ s/DP=//g;
		
		@AFs = ($AFs);
	}
	else {
	
		die ("Unknown caller\n");
	}
	
	return \@Variants, \@AFs, $Depth, $Qual;
}

sub FormatGATKCoordinatesOfVariant {
	
	(my $ReferenceAllele, my $VariantAllele, my $Position) = @_;
	
	while (	substr($ReferenceAllele, -1, 1) eq substr($VariantAllele, -1, 1) ||
			substr($ReferenceAllele, 0, 1) eq substr($VariantAllele, 0, 1)){
	
		if (substr($ReferenceAllele, -1, 1) eq substr($VariantAllele, -1, 1)){
			$ReferenceAllele = substr($ReferenceAllele, 0, (length($ReferenceAllele)-1));
			$VariantAllele = substr($VariantAllele, 0, (length($VariantAllele)-1));
		}
		elsif (substr($ReferenceAllele, 0, 1) eq substr($VariantAllele, 0, 1)){
			$ReferenceAllele = substr($ReferenceAllele, 1, (length($ReferenceAllele)-1));
			$VariantAllele = substr($VariantAllele, 1, (length($VariantAllele)-1));
			$Position++;
		}
	}
	
	if 		($ReferenceAllele eq "")	{$ReferenceAllele 	= "-";}
	elsif	($VariantAllele eq "")		{$VariantAllele 	= "-";}
	
	return ($ReferenceAllele, $VariantAllele, $Position);
}

#
sub ConstructOtherVariant {
					
	(my $Variant, 
	 my $MappingsRef,
	 my $PosToRegionNameRef,
	 my $StrandsRef,
	 my $SDsRef)
	= @_;
	
	(my $ToMapOnChrom, my $ToMapOnStart, my $ToMapOnEnd, my $Change) 	= split(m/~/, $Variant);
	(my $RefNuc, my $AltNuc) 											= split(m/\//, $Change);
	 my $ToMapOnRegionName												= $PosToRegionNameRef->{$ToMapOnChrom}{$ToMapOnStart};
	 my %OtherVariants 													= ();
	 
	#print $Variant . "\t" . $ToMapOnRegionName .  "\t" . $StrandsRef->{$ToMapOnRegionName} .  "\n";
	 
	foreach my $OtherCoord (keys %{$MappingsRef->{"$ToMapOnChrom\~$ToMapOnStart"}}){
		
		(my $OtherChrom, my $OtherStart) 	= split(m/~/, $OtherCoord);
		 my $OtherVariant					= "";
		 my $OtherRegionName				= $PosToRegionNameRef->{$OtherChrom}{$OtherStart};
		 
		 #print \t .  $OtherRegionName . \t . $StrandsRef->{$ToMapOnRegionName} . \n;
		
		### CHECK WHETHER OR NOT THERE IS AN OVERLAPPING SD ###
		### FUNCTION NOT ADOPTED FOR INDELS ###
		 
		if (exists $SDsRef->{$ToMapOnChrom}{$ToMapOnStart}){
		
			foreach my $SD (keys %{$SDsRef->{$ToMapOnChrom}{$ToMapOnStart}}){
			
				(my $SD_Chrom, my $SD_Start, my $SD_End, my $SD_Change) 	= split(m/~/, $SD);
				(my $RefSdNuc, my $AltSdNuc) 								= split(m/\//, $SD_Change);
				
				if ($PosToRegionNameRef->{$SD_Chrom}{$SD_Start} eq $OtherRegionName){
					
					if ($StrandsRef->{$ToMapOnRegionName} eq $StrandsRef->{$OtherRegionName}){
				
						$RefNuc = $AltSdNuc;
					}
					else{
						ReverseComplement($AltSdNuc);
					}
				}
			}	
		}
		
		###########################
		###     SUBSTITUTION    ###
		###########################
		
		if ($RefNuc ne "-" && $AltNuc ne "-"){
			
			my $OtherEnd = $OtherStart;
			
			if ($StrandsRef->{$ToMapOnRegionName} eq $StrandsRef->{$OtherRegionName}){
			
				$OtherVariant = $OtherChrom . "~" . $OtherStart . "~" . $OtherEnd . "~" . $RefNuc . "/" . $AltNuc;
			}
			else {
			
				my $RevComRefNuc = ReverseComplement($RefNuc);
				my $RevComAltNuc = ReverseComplement($AltNuc);
				
				$OtherVariant = $OtherChrom . "~" . $OtherStart . "~" . $OtherEnd . "~" . $RevComRefNuc . "/" . $RevComAltNuc;
			}
		}		
		
		###########################
		### 	 INSERTION 	   ####
		###########################
		
		elsif ($RefNuc eq "-"){
		
			if ($StrandsRef->{$ToMapOnRegionName} eq $StrandsRef->{$OtherRegionName}){
			
				my $OtherEnd = $OtherStart - 1;
			
				$OtherVariant = $OtherChrom . "~" . $OtherStart . "~" . $OtherEnd . "~" . "-/" . $AltNuc;
			}
			else { ### CORRECT ? ###
			
				my $RevComAltNuc = ReverseComplement($AltNuc);
				my $OtherEnd = $OtherStart - 1;
				
				$OtherVariant = $OtherChrom . "~" . $OtherStart . "~" . $OtherEnd . "~" . "-/" . $RevComAltNuc;
			}
		}
		
		###########################
		### 	 DELETION 	   ####
		###########################
		
		else {
		
			if ($StrandsRef->{$ToMapOnRegionName} eq $StrandsRef->{$OtherRegionName}){
				
				my $OtherEnd = $OtherStart + length($RefNuc) - 1;
			
				$OtherVariant = $OtherChrom . "~" . $OtherStart . "~" . $OtherEnd . "~" . $RefNuc . "/" . $AltNuc;
			}
			else {
				
				my $RevComRefNuc 	= ReverseComplement($RefNuc);
				my $RevComAltNuc 	= ReverseComplement($AltNuc);
				my $OtherEnd 		= $OtherStart - length($RefNuc) + 1;
				
				$OtherVariant = $OtherChrom . "~" . $OtherEnd . "~" . $OtherStart . "~" . $RevComRefNuc . "/" . $RevComAltNuc;
			}
		}
		
		$OtherVariants{$OtherVariant} = undef;
	}
		
	return \%OtherVariants;
}

sub ReverseComplement {
    
	(my	$dna) = @_;
     my $revcom = reverse $dna;
		$revcom =~ tr/ACGTacgt/TGCAtgca/;	
     return $revcom;
}