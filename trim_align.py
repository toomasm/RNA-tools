import os
import textwrap
import argparse

from algorithms import perform_trim, perform_alignment

def make_argument_parser():
    '''Returns argument parser for this script.
    '''
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent('''
                                     This script trims sequences from an input file and aligns these to a reference genome.                                     
                                     '''),
                                     fromfile_prefix_chars='@')

    dpg = parser.add_argument_group('Directory and file parameters')
        
    dpg.add_argument('-r', '--reference-genome', 
                     dest='reference_genome', default=None, metavar='S', 
                     help='Reference genome file.')
    
    dpg.add_argument('-i', '--input-file', 
                     dest='input_filename', default=None, metavar='S', 
                     help='Input filename containing fasta sequences.')
        
    return parser

def process_arguments(args):

    # Get the root file name for the ref. genome
    base_name = os.path.basename(args.input_filename)

    # Get the root and extention of the base_name
    root_ext = os.path.splitext(base_name)

    # Make the trimmed base filename based on the read genome name
    trimmed_output_base_filename = root_ext[0] + root_ext[1].replace('.', '_trimmed.')

    # Make alignment filename.
    alignment_output_base_filename = root_ext[0] + '_aligned.sam'

    # Put these in the generation dir...
    trimmed_output_filename = os.path.join('generated', trimmed_output_base_filename)
    alignment_output_filename = os.path.join('generated', alignment_output_base_filename)

    # Perform the trim operation.
    perform_trim(args.input_filename, trimmed_output_filename)

    # Perform the alignment.
    #perform_alignment(args.reference_genome, trimmed_output_filename, alignment_output_filename)

            
if __name__ == '__main__':

    p = make_argument_parser()
    args = p.parse_args()
    process_arguments(args)

