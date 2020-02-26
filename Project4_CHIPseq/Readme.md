# Pipeline to analyse ChIP-seq data 
(From alignment to differential peak calling)

![](./CHIP_Pipeline.jpg)


### Major steps involved in data processing
STEP 1: Download the raw sequencing data (FASTQ files) in local directory

STEP 2: Analyse the quality and statistics of reads using FASTQC-toolkits. 
	[Option: Check for the software installed in cluster/ local machine]

	module load java/1.7.0_25   
	module load fastqc/0.10.1
	fastqc -o /Out_Dir/ -f fastq /Sample.fastq 	# where Sample = fastq samples provided

NOTE: All files provided were of good quality (Phred score > 30)

STEP 3: Align the high quality sequencing reads (from STEP 3) onto human reference genome (hg38) using Hisat (pre installed in Bigred).
This step includes two sub-steps:

	module load hisat/0.1.6        # load module available in cluster/ local machine
	
   1. Building up the indexes for reference genome
   
   	hisat-build -f /Mouse/Mus_musculus.GRCh38.84.dna.toplevel.fa /Mouse/m38.84/m38.84    
	# Refernce genome (and .gtf file for annotation) for mouse was downloaded from Ensembl
		
NOTE: This step is required just once to index the genome build and can be used in next step directly

   2. Alignment:
		
	hisat -p 32 -q -x /Mouse/m38.84/m38.84 -U Rep1.fastq /Output/Rep1.sam

STEP 4: Use samtools (check for module available in cluster/ local machine)
	for post processing of the aligned reads
	
	module load samtools/1.2 						# load module in Bigred
	samtools view -bS /Output/Rep1.sam > /Output/Rep1.bam			# data compression
	samtools cat -o /Output/Mel.bam /Output/Rep1.bam /Output/Rep2.bam	# data concatenation 
	samtools sort /Output/Mel.bam /Output/Mel.sorted	 		# sorting
	samtools index /Output/Mel.sorted.bam					# indexing the bam

STEP 5: Run MACS for chIP seq data analysis and identification of differential binding by following command
	
	module load macs/1.4.2							# load module in Bigred
	cd ./MACS_Output/
	macs14 -t ./ALIGNMENT/CH12.sorted.bam -c ./ALIGNMENT/MEL.sorted.bam -f BAM -g mm -n MELvsCH12 -B -S  
	#Note: -c Control -t Treatment -f format -g organism -n name_study
