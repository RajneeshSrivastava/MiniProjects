## This program is designed to generate script(s) as per the number of fastq files in a given directory.

User can customise the JobTemplate to use this script in any cluster depending on the requirements and modules available

## Current pipeline is following this workflow:

FASTQ ---> SAM ---> BAM ---> SORTED.BAM ---> Transcript_exp.tab + Gene_exp.tab
     Hisat       Samtools           Stringtie

## Built in Unix environment
'''
$python Pipeline_Fastq-to-Quant.py -h
usage: Pipeline_Fastq-to-Quant.py [-h] [-file [fq [fq ...]]]
                                  [-modules [module [module ...]]] [-n N]
                                  [-p P] [-w W] [-index INDEX] [-gtf GTF]
                                  [-odir ODIR]

Routine pipeline for RNA-seq fastq file processing and quant

optional arguments:
  -h, --help            show this help message and exit
  -file [fq [fq ...]]   input fast(q|a) file(s) optional:path/*.fastq
                        (default: None)
  -modules [module [module ...]]
                        module(s) required for data processing (default: None)
  -n N                  number of nodes (default: 1)
  -p P                  number of processors (default: 16)
  -w W                  estimated wall time for process (default: 04:00:00)
  -index INDEX          <path>/reference index file(s) (default:
                        Reference/mouse/index/mm10*)
  -gtf GTF              <path>/reference annotation file (default:
                        Reference/mouse/mm10.gtf)
  -odir ODIR            output directory path (default: )
'''
