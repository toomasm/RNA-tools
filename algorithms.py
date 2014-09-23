import os
import glob

from parse import ParseFastQ

def perform_trim(read_file_path, output_filename, five, three):
    parser = ParseFastQ(read_file_path)    

    print('Performing trim. Making trim file: {}'.format(output_filename))
    output_file = open(output_filename, 'w')

    for record in parser:
        if three == True:
            count = 0
            for x in record[1][::-1]:
                if x == 'A':
                    count += 1
                else:
                    break
            if count != 0:
                record[1] = record[1][:-count]   
                record[3] = record[3][:-count]

        if five == True:
                record[1] = record[1][2:]   
                record[3] = record[3][2:]

        for line in record: 
            output_file.write("%s\n" % line)

def index_genome(reference_genome):
    cmd = "bwa index -a is {0}".format(reference_genome)
    print('Running external command: {}'.format(cmd))
    output = os.system(cmd)

def perform_alignment(reference_genome, index_name, read_filename, output_filename):

    # Test if index files exist.
    index_list = glob.glob('{}*'.format(index_name))
    if len(index_list) < 6:
        print('No index files found. Generating...')
        cmd = "bowtie-build {} {}".format(reference_genome, index_name)
        print('Running external command: {}'.format(cmd))
        output = os.system(cmd)

    # Run the alignment.
    cmd = "bowtie -v 1 -k 7 -S -p 2 --best --strata {} {} {}".format(index_name, read_filename, output_filename)
    print('Running external command: {}'.format(cmd))
    output = os.system(cmd)

    print('Created alignment file: {}'.format(output_filename))

def sam_to_bam(sam_file, bam_file):
    cmd = "samtools view -Sb  {}  >  {}".format(sam_file, bam_file)
    print('Running external command: {}'.format(cmd))
    output = os.system(cmd)

def check_file_type(file_extension):
    if file_extension == '.fastq' or file_extension == '.fastq.bz2' or file_extension == '.fastq.gz':
        type = 'fastq'
    elif file_extension == '.sam':
        type = 'sam'
    elif file_extension == '.bam':
        type = 'bam'
    else:
        type = 'unknown'
    return type