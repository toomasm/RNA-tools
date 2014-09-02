import os
import textwrap
import argparse

from trimmer import five_prime_trimmer, three_prime_trimmer, three_prime_checker
from align_and_index import perform_alignment, sam_to_bam
from index_generator import generate_index
def make_argument_parser():
    '''Returns argument parser for this script.
    '''
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                             description=textwrap.dedent('''
                             This script takes a fastq file as an input, trims it as required and alignes it to reference genome using bowtie.
                             Finally it generates an index file containing read names as indexes. The columns represent different ribosomal RNA operons
                             in E. coli. This script adds a mapped reads 5' or 3' position to the appropiate rRNA operon (one read can be aligned to
                             multiple locations). Also the orienation of the strand is taken into the account.                                      
                             '''))

    trimmer_args = parser.add_argument_group('fastq file trimming parameters')

    trimmer_args.add_argument('-5','--N-from-5p-to-trim',
                              default=None, type=int, dest='N_from_5p_to_trim',
                              help='Trims your specifified number (takes INT as argument) of nucleotides from the five prime end of every read')

    trimmer_args.add_argument('-3','--trim-all-A-from-3p',
                              default=None, dest = 'trim_all_A_from_3p', action='store_true',
                              help='If True trims all "A" nucleotides from the 3 prime end of the fastq reads.\
                              It removes "A" nucleotides till it reaches a none A nucleotide ')

    trimmer_args.add_argument('-c','--three-prime-check',
                              default=None, dest ='three_prime_check', action='store_true',
                              help="If True divides the reads of a sam file between two new files according to nucleotide following\
                              the mapped three prime end of the read. If the nucleotide is not A the read is true and is saved in the\
                              correct reads file. If the nucleotide is A we have no means to know if the read has been correctly trimmed or not\
                              and it is saved in shady reads file.")

    data_file_args = parser.add_argument_group('Data files descriptions and parameters')
    
    data_file_args.add_argument('-r', '--reference-genome', 
                                dest='reference_genome', default=None, 
                                help='Input a reference genome file with .fasta extension. ')
    
    data_file_args.add_argument('-i', '--input-file', 
                                dest='input_filename', default=None, 
                                help='Input fastq filename')

    aligment_args = parser.add_argument_group('Aligment and indexing arguments')

    aligment_args.add_argument('-a', '--alignment',
                                default=None, dest='alignment', action="store_true",
                                help='Alignes your reads using bowtie against reference genome.\
                                Converts the sam to bam and deletes the original sam file.\
                                output is SAM file' )

    aligment_args.add_argument('--index', action='store_true',
                               dest='index_the_fasta',
                               help='Indexes reference genome with bowtie indexer' )

    aligment_args.add_argument('-b', '--sam-to-bam', action ='store_true',
                               dest=('sam_to_bam'), default = None,
                               help='Converts sam to bam')
        
    return parser

def process_arguments(args):

    # Get the root file name for the input file and reference genome fasta file.
    base_name = os.path.basename(args.input_filename)

    # Get the root and extention of the base_name and the base_name_ref_gen.
    root_ext = os.path.splitext(base_name)
   
    # Make the trimmed base filename based on the read genome name.
    trimmed_output_base_filename = root_ext[0] + root_ext[1].replace('.', '_trimmed.')
    
    # Make alignment filename for sam and bam files.
    SAM_output_base_filename = root_ext[0] + '_aligned.sam'
    BAM_output_base_filename = root_ext[0] + '_aligned.bam'
    index_file_base_name = root_ext[0] + '.csv'

    # Put trimmed files in the generated dir, and alignement files into aligned dir. If such dir not present, create one.
    if not os.path.exists('aligned'):
        os.makedirs('aligned')
        
    if not os.path.exists('generated'):
        os.makedirs('generated')
        
    trimmed_output_filename = os.path.join('generated', trimmed_output_base_filename)
    SAM_output_filename = os.path.join('aligned', SAM_output_base_filename)
    BAM_output_filename = os.path.join('aligned', BAM_output_base_filename)
    index_filename = os.path.join('aligned', index_file_base_name)
    #Make genome index (bowtie) filename
    if not args.reference_genome == None:
        base_name_ref_gen = os.path.basename(args.reference_genome)
        root_ext_ref_gen = os.path.splitext(base_name_ref_gen)
        reference_genome_index = root_ext_ref_gen[0]
    # index_list = glob.glob('{}*'.format(reference_genome_index_base)).
    
    
    # Perform the trim operation.
    if root_ext[-1] == '.fastq':
        if args.N_from_5p_to_trim >= 1  or args.trim_all_A_from_3p == True:
            five_prime_trimmer(args.input_filename, trimmed_output_filename, args.N_from_5p_to_trim, args.trim_all_A_from_3p)
        else:            
            print "No trimming"
    
    if args.alignment == True and root_ext[-1] == '.fastq':
        perform_alignment(args.reference_genome, reference_genome_index, args.input_filename, SAM_output_filename)
        sam_to_bam(SAM_output_base_filename, BAM_output_filename)
        generate_index()
        
    elif root_ext[-1] != '.fastq' and root_ext_ref_gen[-1] != '.fasta':
        print "Your input file has to be a fastq file (-i FILENAME.fasta) and reference genome file a fasta file (-r FILENAME.fasta)"
       
    if args.sam_to_bam == True and root_ext[-1]== '.sam':
        sam_filename_path = "aligned/" + SAM_output_base_filename
        sam_to_bam(sam_filename_path, BAM_output_filename)
        bam_filename_path = "aligned/" + BAM_output_base_filename
        generate_index(bam_filename_path, index_filename)
        
        
    # Perform the alignment.
    #perform_alignment(args.reference_genome, trimmed_output_filename, alignment_output_filename).


if __name__ == '__main__':

    p = make_argument_parser()
    args = p.parse_args()
    process_arguments(args)

