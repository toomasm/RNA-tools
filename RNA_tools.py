import os
import textwrap
import argparse


def make_argument_parser():
    '''Returns argument parser for this script.
    '''
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent('''
                                     This script unifies our tools for trimming, modifying and aligning our reads and 
                                     also analyzes mapped reads and generates an index file for plotting the data.                                     
                                     '''),
                                     fromfile_prefix_chars='@')

    dpg = parser.add_argument_group('Directory and file parameters')
        
    dpg.add_argument('-i1', '--input-file', 
                     dest='input_one', default=None, 
                     help='Can be either a .fastq, .bam or .csv file.')
    
    dpg.add_argument('-i2', '--second-input-file', 
                     dest='input_two', default=None, 
                     help='optional .fasta file')

    dpg.add_argument('-3', '--trim-3\'-end', 
                     dest='trim_3', default=None, action = 'store_true', 
                     help='If True, program will trim 2 nucleotides from 3\' end.')

    dpg.add_argument('-5', '--trim-5\'-end', 
                     dest='trim_5', default=None, action = 'store_true', 
                     help='if True, program will trim A-s from 5\' end.')

    dpg.add_argument('-a', '--align reads to reference genome', action = 'store_true',
                     dest='align', default=None, 
                     help='If True, program will align reads to reference genome.')

    dpg.add_argument('-i', '--generate index file', action = 'store_true',
                     dest='index', default=None, 
                     help='If True, program will generate index file according to which operons a certain read aligns to.')

    return parser

import algorithms
import index_generator

def process_arguments(args):

    # Get the root file name for first input file
    base_name_one = os.path.basename(args.input_one)

    # Get the root and extention of base_name_one
    root_ext_one = os.path.splitext(base_name_one)

    # Make the trimmed reads file output name
    trimmed_output_filename = os.path.join('generated', (root_ext_one[0] + root_ext_one[1].replace('.', '_trimmed.')))

    # Make aligned output name
    aligned_sam_name = os.path.join('generated', (root_ext_one[0] + '_aligned.sam'))
    aligned_bam_name = aligned_sam_name.replace('.sam', '.bam')

    # Get the root file name for second input file
    base_name_two = os.path.basename(args.input_two)

    # Get the root and extension of base_name_two
    root_ext_two = os.path.splitext(base_name_two)

    # Make genome index name
    genome_index_name = 'generated/' + root_ext_two[0]

    # Make index name
    index_name = 'generated/' + (root_ext_one[0] + '_index.csv')

    # Created a directory for working files
    if not os.path.exists('generated'):
        os.makedirs('generated')

    # Perform trim
    if root_ext_one[1] == '.fastq':
        if args.trim_3 or args.trim_5:
            algorithms.perform_trim(args.input_one, trimmed_output_filename, args.trim_5, args.trim_3)
    elif root_ext_one[1] != '.fastq' and (args.trim_3 or args.trim_5):
        print 'Input file is not a .fastq file'

    # Perform alignment
    if args.align:
        if args.trim_3 or args.trim_5:
            algorithms.perform_alignment(args.input_two, genome_index_name, trimmed_output_filename, aligned_sam_name)
        else:
            algorithms.perform_alignment(args.input_two, genome_index_name, args.input_one, aligned_sam_name)
        algorithms.sam_to_bam(aligned_sam_name, aligned_bam_name)    

    # Generate index file
    if args.index:
        if args.align:
            index_generator.generate_index(aligned_bam_name, index_name)
        else:
            index_generator.generate_index(args.input_one, index_name)



if __name__ == '__main__':

    p = make_argument_parser()
    args = p.parse_args()
    process_arguments(args)