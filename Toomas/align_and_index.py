import os
import textwrap
import glob

def perform_alignment(reference_genome, reference_genome_index_base, fastq_filename, sam_output_filename):
    #bowtie -v 1 -k 7 -p 2 -S --trim5 2 --best --strata MG1655_913.3 MazF2h_PNK.fastq MazF_PNK.sam
    # Test if index files exist. Creates a index file for reference genome if it does not.
    index_list = glob.glob('{}*'.format(reference_genome_index_base))
    if len(index_list) != 7:
        #Creates a index file for reference genome if it does not.
        print('No index files found. Generating...')
        cmd = "bowtie-build {0} {1} ".format(reference_genome, reference_genome_index_base)
        print('Running external command: {0}'.format(cmd))
        output = os.system(cmd)

    # Map the reads to reference genome using bowtie 1.
    cmd = "bowtie -v 1 -k 7 -p 2 -S --best --strata {0} {1} {2} ".format(reference_genome_index_base, fastq_filename_path, sam_output_filename)
    print('Running external command: {}'.format(cmd))
    output = os.system(cmd)

    print('Created sam alignment file: {}'.format(sam_output_filename))

def sam_to_bam(sam_file, bam_file):
    #Convert sam to bam      
    cmd = 'samtools view -b -S {0} > {1}'.format(sam_file, bam_file)
    output = os.system(cmd)

    print('Created bam alignment file: {}'.format(bam_file))