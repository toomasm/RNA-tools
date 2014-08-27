import os
import textwrap
import argparse

import pysam
from Bio import SeqIO

def make_argument_parser():
    '''Returns argument parser for this script.
    '''
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent('''
                                     This script sorts aligned reads according to whether their 3' end has been trimmed correctly or there is cause for suspicion.                                     
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

    # Make the correctly aligned reads file output name
    correct_output_base_filename = root_ext[0] + root_ext[1].replace('.', '_correct.')

    # Make the suspiciously aligned reads file output name
    suspicious_output_base_filename = root_ext[0] + root_ext[1].replace('.', '_suspicious.')

    # Put these in the generation dir...
    correct_output_filename = os.path.join('generated_3prim', correct_output_base_filename)
    suspicious_output_filename = os.path.join('generated_3prim', suspicious_output_base_filename)

    # Perform the sorting operation.
    sort_samfile(args.input_filename, args.reference_genome, correct_output_filename, suspicious_output_filename)


def sort_samfile(samfile, reference_genome, correct_output_filename, suspicious_output_filename):
  
    for seq_record in SeqIO.parse(reference_genome, "fasta"):
        sequence=seq_record.seq

    samfile = pysam.Samfile(samfile, "r" )
    header_line = samfile.header
    text_line = samfile.text

    output_correct=pysam.Samfile(correct_output_filename, 'w', header = header_line, text = text_line)
    output_suspicious=pysam.Samfile(suspicious_output_filename, 'w', header = header_line, text = text_line)

    for read in samfile:
        if read.is_reverse == True:
            if read.pos == sequence[0]:
                check_pos = sequence[-1]
            else:
                check_pos = read.pos - 1
                if sequence[check_pos] == 'A':
                    output_suspicious.write(read)
                else:
                    output_correct.write(read) 
        else:
            check_pos = read.pos + read.qlen
            if check_pos == len(sequence):
                if sequence[0] == 'A':
                    output_suspicious.write(read)
                else:
                    output_correct.write(read)
            elif sequence[check_pos] == 'A':
                output_suspicious.write(read)
            else:
                output_correct.write(read)
            

if __name__ == '__main__':

    p = make_argument_parser()
    args = p.parse_args()
    process_arguments(args)