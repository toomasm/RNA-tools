import os
import textwrap
import argparse

from trimmer import read_trimmer, three_prime_checker
from align_and_index import perform_alignment, sam_to_bam
from index_generator_2 import generate_index

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
                              default=0, type=int, dest='N_from_5p_to_trim',
                              help='Trims your specifified number (takes INT as argument) of nucleotides from the five prime end of every read')

    trimmer_args.add_argument('-3','--trim-all-A-from-3p',
                              default=False, dest = 'trim_all_A_from_3p', action='store_true',
                              help='If True trims all "A" nucleotides from the 3 prime end of the fastq reads.\
                              It removes "A" nucleotides till it reaches a none A nucleotide.\
                              \
                              Also if set, tells index generator that it is dealing with three prime mapped reads.')

    trimmer_args.add_argument('-c','--three-prime-check',
                              default=False, dest ='three_prime_check', action='store_true',
                              help="If True divides the reads of a bam file between two new files according to nucleotide following\
                              the mapped three prime end of the read. If the nucleotide is not A the read is true and is saved in the\
                              correct reads file. If the nucleotide is A we have no means to know if the read has been correctly trimmed or not\
                              and it is saved in shady reads file.")

    data_file_args = parser.add_argument_group('Data files descriptions and parameters')
    
    data_file_args.add_argument('-r', '--reference-genome', 
                                dest='reference_genome', default=None, 
                                help='Input a reference genome file with .fasta extension. ')
    
    data_file_args.add_argument('-i', '--input-file', 
                                dest='input_filename', default=None, required=True,
                                help='Input filename. Can be fastq, SAM or BAM format.\
                                Required argument')

    alignment_args = parser.add_argument_group('Aligment and indexing arguments')

    alignment_args.add_argument('-a', '--alignment',
                                default=None, dest='alignment', action="store_true",
                                help='Alignes your reads using bowtie against reference genome.\
                                Converts the sam to bam and deletes the original sam file.\
                                output is SAM file' )

    alignment_args.add_argument('--index', action='store_true',
                               dest='index_the_bam',
                               help='Indexes the bam files reads. Checks to which rRNA operons (E. coli MG1655 v 913.3) it maps to.')
        
    return parser

def trimming_operations(input_filename, root_ext, N_from_5p_to_trim, trim_all_A_from_3p):
    '''Function that checks if trimming options have been set and conditions met. If yes, the reads are trimmed accordingly.
       Returns trimmed file path. If not trimmed, returns input fastq file path.'''

    #Checks if your 5 prime trimming option entry is valid.
    
    if type(N_from_5p_to_trim) != int:
        raise UserWarning ('Input has to be a single number. This is not a single number!')
    
    #Check if args for trimming have been set.
    if N_from_5p_to_trim >= 1  or trim_all_A_from_3p == True:
        # Make filename for trimmed fastq file. Put it in "trimmed"" directory.
        trimmed_output_filename = os.path.join('trimmed', (root_ext[0] + root_ext[-1].replace('.', '_trimmed.')))
        #Perform the trim operation.
        read_trimmer(input_filename, trimmed_output_filename, N_from_5p_to_trim, trim_all_A_from_3p)
        return trimmed_output_filename
    elif  N_from_5p_to_trim == 0  and trim_all_A_from_3p == False:
        print "No trimming of the reads"
        return args.input_filename

def alignment_operations(input_filetype, reference_genome, root_ext, fastq_filename):
    '''Processes alignment operations. Checks if input and reference options and filetypes are valid.
       If yes performs alignment and returns path to alignment SAM file. '''

    #Checks wether reference genome option is set.
    if reference_genome == None:
        raise UserWarning('You need to set reference_genome option (-r).'\
                          'It takes ".fasta" file as value')

    #Make genome index (bowtie) filename
    #Get root filename for reference genome.
    base_name_ref_gen = os.path.basename(reference_genome)
    #Get root and extention of the reference genome fasta file.
    root_ext_ref_gen = os.path.splitext(base_name_ref_gen)
    #Ger the base of reference genome name.
    reference_genome_index = root_ext_ref_gen[0]
    #Get the file extention of reference genome.
    reference_genome_filetype = root_ext_ref_gen[-1]

    #Checks wether input and reference filetypes are valid. If not, UserWarning is raised.
    if input_filetype != '.fastq' or reference_genome_filetype != '.fasta':
        raise UserWarning('Your input file has to be a fastq file (-i FILENAME.fasta) ' \
                          'and reference genome file a fasta file (-r FILENAME.fasta).')

    # Makes alignment filename for SAM file. The file will be put in "aligned" directory.
    SAM_output_filename = os.path.join('aligned', (root_ext[0] + '_aligned.sam'))
    print "Generating SAM file using Bowtie."
    #Generates a SAM file using Bowtie algorithm.
    perform_alignment(reference_genome, reference_genome_index, fastq_filename, SAM_output_filename)
    return SAM_output_filename

def BAM_maker_operations(SAM_output_filename):
    '''Converts SAM to BAM using SAMtools.'''

    #Makes name for BAM output file.
    BAM_output_filename = SAM_output_filename.replace('.sam', '.bam')
    #Converts SAM to BAM using Samtools.
    print "Converting SAM to BAM"
    sam_to_bam(SAM_output_filename, BAM_output_filename)
    return BAM_output_filename

def dividing_3prime_seq_aligned_reads(BAM_output_filename, reference_genome, bam_list, root_ext):
    '''Function that makes necessary operations for dividing 3' polyA trimming for 3' end seq reads.
       Divides reads aligned in BAM file to ok and shady reads files. This is used only in case 3' polyA
       trimming for 3' seq reads. It checks whether whether after trimming the nucleotide following the
       mapped read is A or not. When it is A it is no way to be sure if the mapped reads 3' position is
       correct and read is added to shady reads bam file. Otherwise read is added to ok reads file.'''
    
    print os.getcwd() + '/aligned'
    print os.path.split(os.path.abspath(BAM_output_filename))[0]
    print os.path.join(os.getcwd(), 'aligned')
    #Checks if the BAM files path is the same as our created aligned BAM files path.
    if os.path.split(os.path.abspath(BAM_output_filename))[0] == os.path.join(os.getcwd(), 'aligned'):
        print 'Fuck Yeah!'
    else:
        print 'Nope!'
    
    BAM_base = os.path.basename(BAM_output_filename)
    print BAM_base
    
    #Checks if BAM file is named according to our codes rules.
    #Makes a new name for ok and shady reads BAM files accordingly.
    if '_aligned.bam' in BAM_base:
        BAM_output_filename_ok = os.path.join('aligned', ('ok_' + BAM_base))
        BAM_output_filename_shady = os.path.join('aligned', ('shady_' + BAM_base))
    else:
        BAM_output_filename_ok = os.path.join('aligned', ('ok_' + root_ext[0] + '_aligned.bam'))
        BAM_output_filename_shady = os.path.join('aligned', ('shady_' + root_ext[0] + '_aligned.bam'))
        
    #Divides three prime end trimmed reads into ok and suspicious reads.
    three_prime_checker(BAM_output_filename, reference_genome, BAM_output_filename_ok, BAM_output_filename_shady)
    #Adds these bam files to bam_list for further processing.
    bam_list.extend([BAM_output_filename_ok, BAM_output_filename_shady])
    return bam_list
    
def index_operations(root_ext, bam_list, base_name, three_prime):
    # Makes index files from BAM files.

    for BAM in bam_list:
        #Generates index file name if BAM file names were created by or named similarly to current program.
        if '_aligned.bam' in BAM:
            index_filename = BAM.replace('_aligned.bam', '_index.csv')
        #If BAM is named differently then just replaces '.bam' with '_index.csv'.
        else:
            index_filename = BAM.replace('.bam', '_index.csv')

        # Makes an index file from BAM. If three_prime check is True makes the indexer count three prime ends instead of five prime.      
        generate_index(BAM, index_filename, three_prime)
        print 'Generated index named {0}'.format(os.path.basename(index_filename))
        
def process_arguments(args):

    abspath = os.path.abspath(args.input_filename)
    dir_name = os.path.dirname(abspath)
    os.chdir(dir_name)
    
    # Get the root file name for the input file.
    base_name = os.path.basename(args.input_filename)

    # Get the root and extention of the base_name.
    root_ext = os.path.splitext(base_name)

    #Get the extention of input file.
    input_filetype = root_ext[-1]

    #Sets the value of SAM_output_filename as empty string.
    #Acts as a flag and if SAM-s present stores SAM file path.
    SAM_output_filename = ''

    #Creates a list for BAM files. This list will contain BAM files to be further processed.
    bam_list = []

    #If not present, generate dir for trimmed files and dir for aligned and index files.
    if not os.path.exists('aligned'):
        os.makedirs('aligned')
        
    if not os.path.exists('trimmed'):
        os.makedirs('trimmed')

    #Checks if input file is '.fastq'.
    if input_filetype == '.fastq':
        #Function that checks if trimming options have been set and conditions met. If yes, the reads are trimmed accordingly.
        #Returns trimmed file path. If not trimmed, returns input fastq file path.
        fastq_filename = trimming_operations(args.input_filename, root_ext,
                                                      args.N_from_5p_to_trim, args.trim_all_A_from_3p)
        #Checks if alignement option has been set.
        if args.alignment == True:
            #If alignment option has been set, checks whether conditions (filetypes, options etc.) have been met.
            #If all is ok, performs alignment.Output is path to alignment file.
            SAM_output_filename = alignment_operations(input_filetype,args.reference_genome,
                                                       root_ext, fastq_filename)
    
    #If SAM is the input file assignes its path to SAM_output_filename variable.
    if input_filetype == '.sam':
        SAM_output_filename = input_filename
        #This flag will be used for determining what operations to skip.

    if SAM_output_filename:
        #Converts SAM to BAM. Returns BAM output path.
        BAM_output_filename = BAM_maker_operations(SAM_output_filename)
        #Checks if the reads were 3' polyA trimmed.
        #If yes, checks for possible false trimming and divides them to ok and shady reads.
        if args.three_prime_check == True:
            #Output is a list of bamfiles for further processing.
            bam_list = dividing_3prime_seq_aligned_reads(BAM_output_filename, args.reference_genome, bam_list, root_ext)
        else:
            bam_list.append(BAM_output_filename)

    #If input file is BAM, adds it to bam_list for further processing.
    if input_filetype == '.bam':
        if args.three_prime_check == True:
            #Output is a list of bamfiles for further processing.
            bam_list = dividing_3prime_seq_aligned_reads(args.input_filename, args.reference_genome, bam_list, root_ext)
        else:
            bam_list.append(args.input_filename)

    #Makes index files from BAM files in bam_list if argument for generating the index is set True. 
    if args.index_the_bam == True and bam_list:
        index_operations(root_ext, bam_list, base_name, args.trim_all_A_from_3p)
   
if __name__ == '__main__':

    p = make_argument_parser()
    args = p.parse_args()
    process_arguments(args)

