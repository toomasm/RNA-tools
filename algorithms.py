import os
import glob

from parse import ParseFastQ

def perform_trim(read_file_path, output_filename):
    parser = ParseFastQ(read_file_path)    

    print('Performing trim. Making trim file: {}'.format(output_filename))
    output_file = open(output_filename, 'w')

    for record in parser:
        count = 0
        for x in record[1][::-1]:
            if x == 'A':
                count += 1
            else:
                break
        if count == 0:
            record[1] = record[1][2:]   
            record[3] = record[3][2:]
        else: 
            record[1] = record[1][2:-count]   
            record[3] = record[3][2:-count]

        for line in record: 
            output_file.write("%s\n" % line)

def index_genome(reference_genome):
    cmd = "bwa index -a is {0}".format(reference_genome)
    print('Running external command: {}'.format(cmd))
    output = os.system(cmd)

def perform_alignment(reference_genome, read_filename, output_filename):

    # Test if index files exist.
    index_list = glob.glob('{}*'.format(reference_genome))
    if len(index_list) != 6:
        print('No index files found. Generating...')
        cmd = "bwa index {}".format(reference_genome)
        print('Running external command: {}'.format(cmd))
        output = os.system(cmd)

    # Run the alignment.
    cmd = "bwa mem {} {} > {}".format(reference_genome, read_filename, output_filename)
    print('Running external command: {}'.format(cmd))
    output = os.system(cmd)

    print('Created alignment file: {}'.format(output_filename))
