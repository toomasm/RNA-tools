import os
import subprocess
import pandas as pd
import pysam
from sets import Set
from progressbar import AnimatedMarker, Bar, BouncingBar, Counter, ETA, \
    ProgressBar, ReverseBar, RotatingMarker, \
    SimpleProgress, Timer


def generate_index(bam_input_filename):

    number_of_reads = int(subprocess.check_output(["samtools", "view", "-c", bam_input_filename]))

    widgets = [Bar('>'), ' ', ETA(), ' ', ReverseBar('<')]
    pbar = ProgressBar(widgets=widgets, maxval=number_of_reads).start()
    i = 0  

    samfile = pysam.Samfile(bam_input_filename, 'rb')
    names = set()
    columns = ['rrsA','rrsB','rrsC','rrsD','rrsE','rrsG','rrsH','P1','P2','P3','P4','P5','P6','P7']

    for read in samfile:
        names.add(read.qname)
        i += 1
        pbar.update(i)
    pbar.finish()
         
    names_list = list(names)

    df = pd.DataFrame(index=names_list, columns=columns)


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
                if (position > 4035153 and position < 4040906):
                    df['rrsA'][read.qname] = 1
                    df['P1'][read.qname] = position_marker
                elif (position > 4165428 and position < 4172057):
                    df['rrsB'][read.qname] = 1
                    df['P2'][read.qname] = position_marker
                elif (position > 3941327 and position < 3946872):
                    df['rrsC'][read.qname] = 1
                    df['P3'][read.qname] = position_marker
                elif (position > 3423194 and position < 3429236):
                    df['rrsD'][read.qname] = 1
                    df['P4'][read.qname] = position_marker
                elif (position > 4207532 and position < 4213234):
                    df['rrsE'][read.qname] = 1
                    df['P5'][read.qname] = position_marker
                elif (position > 2725746 and position < 2731600):
                    df['rrsG'][read.qname] = 1
                    df['P6'][read.qname] = position_marker
                elif (position > 223408 and position < 229167):
                    df['rrsH'][read.qname] = 1
                    df['P7'][read.qname] = position_marker
            elif read.is_reverse == True:
                position = (read.pos + read.rlen)
                position_marker = '-' + str(read.pos + read.rlen)
                if (position > 4035153 and position < 4040906):
                    df['rrsA'][read.qname] = 1
                    df['P1'][read.qname] = position_marker
                elif (position > 4165428 and position < 4172057):
                    df['rrsB'][read.qname] = 1
                    df['P2'][read.qname] = position_marker
                elif (position > 3941327 and position < 3946872):
                    df['rrsC'][read.qname] = 1
                    df['P3'][read.qname] = position_marker
                elif (position > 3423194 and position < 3429236):
                    df['rrsD'][read.qname] = 1
                    df['P4'][read.qname] = position_marker
                elif (position > 4207532 and position < 4213234):
                    df['rrsE'][read.qname] = 1
                    df['P5'][read.qname] = position_marker
                elif (position > 2725746 and position < 2731600):
                    df['rrsG'][read.qname] = 1
                    df['P6'][read.qname] = position_marker
                elif (position > 223408 and position < 229167):
                    df['rrsH'][read.qname] = 1
                    df['P7'][read.qname] = position_marker
                
    pbar.finish()        

    samfile.close()
    df.to_csv('index.csv',index=True,header=True)

if __name__ == '__main__':
    bam = '/media/markus/DATA2/MazF_PNK_sorted.bam'
    generate_index(bam)