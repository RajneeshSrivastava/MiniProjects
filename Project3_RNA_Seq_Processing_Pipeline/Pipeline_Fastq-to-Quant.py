import os
import argparse


#Running a Job on HPC using PBS. Flexible with any cluster. Just change the JobTemplate!
JobTemplate = """
#!/bin/bash
#PBS -l nodes=%d:ppn=%d 
#PBS -l walltime=%s 
#PBS -l gres=ccm
#PBS -N %s"""

'''
This program is designed to generate script(s) as per the number of fastq files in a given directory.
User can customise the JobTemplate to use this script in any cluster depending on the requirements and modules available
Current pipeline is following this workflow:
FASTQ ---> SAM ---> BAM ---> SORTED.BAM ---> Transcript_exp.tab + Gene_exp.tab
     Hisat       Samtools           Stringtie
'''

def main():
    parser = argparse.ArgumentParser(
                    description='Routine pipeline for RNA-seq fastq file processing and quant',
                    formatter_class=argparse.ArgumentDefaultsHelpFormatter
                                    )
    parser.add_argument('-file', metavar='fq', type = str, nargs='*', help='input fast(q|a) file(s) optional:path/*.fastq')
    parser.add_argument('-modules', metavar='module', nargs='*', help='module(s) required for data processing')
    parser.add_argument('-n', type = int, default = 1, help='number of nodes')
    parser.add_argument('-p', type = int, default = 16, help='number of processors')
    parser.add_argument('-w', type = str, default = '04:00:00', help='estimated wall time for process')
    parser.add_argument('-index', type = str, default = 'Reference/mouse/index/mm10*', help='<path>/reference index file(s)')
    parser.add_argument('-gtf', type = str, default = 'Reference/mouse/mm10.gtf', help='<path>/reference annotation file')
    parser.add_argument('-odir', type = str, default = '', help='output directory path')
    args = parser.parse_args()
    
    for x in args.file:
        x_name = x.split('/')[-1].split('.')[0]
        x_dir = os.path.dirname(x)
        with open(x + '.sh', 'w') as f:
            f.write(JobTemplate %(args.n,args.p,args.w,x_name)+'\n\n')     #os.path.dirname(x) FOR fastq directory 
            for m in args.modules:
                f.write('module load '+m+'\n')
            f.write('\ncd '+x_dir+'\n\n')
            '''
            ALIGNER
            User can add multiple aligner here and make use of it by tweaking in the below code as per the user defined paramenters in tool guidelines:
            if 'aligner' in str(args.modules):
                f.write('ccmrun aligner -p '+str(args.p)+ ' -q -x '+str(args.index)+ ' -u '+x+ ' -S ' + args.odir+x_name+'.sam'+'\n')
            '''
            if 'hisat' in str(args.modules):
                f.write('ccmrun hisat -p '+str(args.p)+ ' -q -x '+str(args.index)+ ' -u '+x+ ' -S ' + args.odir+x_name+'.sam'+'\n')
            '''
            SAMTools is used fpr post processing of aligned reads. This tool can be customised as per user choice.
            In case, the aligner is providing the processed reads, samtool module can be ignored. Example - TopHat!            
            '''      
            if 'samtools' in str(args.modules):
                f.write('ccmrun samtools view -bS '+args.odir+x_name+'.sam > '+ args.odir+x_name+'.bam'+'\n')
                f.write('ccmrun samtools sort '+args.odir+x_name+'.bam '+ args.odir+x_name+'.sorted'+'\n')
                f.write('ccmrun samtools index '+args.odir+x_name+'.sorted.bam'+'\n')
            '''
            QUANTIFIER
            User can add multiple quantification tool and make use of it by tweaking in the below code as per the user defined paramenters in tool guidelines:
            '''
            if 'stringtie' in str(args.modules):
                f.write('ccmrun stringtie '+args.odir+x_name+'.sorted.bam -p '+str(args.p)+ ' -G '+str(args.gtf)+' -e -o '+args.odir+x_name+'.gtf'+' -A '+args.odir+x_name+'.tab'+'\n')

    return
if __name__ == "__main__":
    main()
