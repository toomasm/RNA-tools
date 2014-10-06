import os
import subprocess
import pandas as pd
import pysam
from collections import namedtuple, OrderedDict
from sets import Set
from progressbar import AnimatedMarker, Bar, BouncingBar, Counter, ETA, \
    ProgressBar, ReverseBar, RotatingMarker, \
    SimpleProgress, Timer


def generate_index(bam_input_filename, index_filename):

    number_of_reads = int(subprocess.check_output(["samtools", "view", "-c", bam_input_filename]))

    widgets = [Bar('>'), ' ', ETA(), ' ', ReverseBar('<')]
    pbar = ProgressBar(widgets=widgets, maxval=number_of_reads).start()
    i = 0  

    GeneInfo = namedtuple('GeneInfo', ['start_pos', 'end_pos', 'column_label'])

    gene_dic = OrderedDict([('rrsA', GeneInfo(4035153, 4040906, 'P1')), 
                            ('rrsB', GeneInfo(4165428, 4172057, 'P2')),
                            ('rrsC', GeneInfo(3941327, 3946872, 'P3')),
                            ('rrsD', GeneInfo(3423194, 3429236, 'P4')),
                            ('rrsE', GeneInfo(4207532, 4213234, 'P5')),
                            ('rrsG', GeneInfo(2725746, 2731600, 'P6')),
                            ('rrsH', GeneInfo(223408, 229167, 'P7'))])

    plabels = [d.column_label for g, d in gene_dic.items()]
    columns = ['Readnames'] + gene_dic.keys() + plabels

    samfile = pysam.Samfile(bam_input_filename, 'rb')
    names = set()
    for read in samfile:
        names.add(read.qname)
        i += 1
        pbar.update(i)
    pbar.finish()
         
    names_list = list(names)

    df = pd.DataFrame(columns=columns)

    widgets = [Bar('>'), ' ', ETA(), ' ', ReverseBar('<')]
    pbar = ProgressBar(widgets=widgets, maxval=number_of_reads).start()
    i = 0 

    samfile = pysam.Samfile(bam_input_filename, 'rb')

    for read in samfile:
        i += 1
        pbar.update(i)
        if read.is_unmapped == False:
            if read.is_reverse == False:
                position = (read.pos+1)
                position_marker = (read.pos+1)
                for gene, gene_data in gene_dic.items():
                    if (position > gene_data.start_pos and position < gene_data.end_pos):
                        df[gene][read.qname] = 1
                        df[gene_data.column_label][read.qname] = position_marker
            elif read.is_reverse == True:
                position = (read.pos + read.rlen)
                position_marker = '-' + str(read.pos + read.rlen)
                for gene, gene_data in gene_dic.items():
                    if (position > gene_data.start_pos and position < gene_data.end_pos):
                        df[gene][read.qname] = 1
                        df[gene_data.column_label][read.qname] = position_marker
    pbar.finish()        

    samfile.close()
    df.to_csv(index_filename,index=False,header=True)

