import matplotlib
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
import pandas as pd
from collections import namedtuple, OrderedDict, Counter, defaultdict
import math
import numpy as np
from Bio import SeqIO

GeneInfo = namedtuple('GeneInfo', ['start_pos', 'end_pos', 'column_label'])

gene_dic = OrderedDict([('rrsA', GeneInfo(4035153, 4040906, 'P1')),
                        ('rrsB', GeneInfo(4165428, 4172057, 'P2')),
                        ('rrsC', GeneInfo(3941327, 3946872, 'P3')),
                        ('rrsD', GeneInfo(3423194, 3429236, 'P4')),
                        ('rrsE', GeneInfo(4207532, 4213234, 'P5')),
                        ('rrsG', GeneInfo(2725746, 2731600, 'P6')),
                        ('rrsH', GeneInfo(223408, 229167, 'P7'))])

all_operons = ['rrsA', 'rrsB', 'rrsC', 'rrsD', 'rrsE', 'rrsG', 'rrsH']
operon_positions = ['P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7']


#df1 = '/home/toomas/index_test.csv'
df1 = '/media/toomas/DATA1/aligment_bowtie/MazF_2h/aligned/MazF2h_PNK_index.csv'
df2 = '/media/toomas/DATA1/aligment_bowtie/MazF_2h/aligned/MazF2h_none_index.csv'
#df1 = pd.read_csv('/media/toomas/DATA1/aligment_bowtie/MazF_2h/aligned/MazF2h_PNK_index.csv)
#df2 = pd.read_csv('/media/toomas/DATA1/aligment_bowtie/MazF_2h/aligned/MazF2h_none_index.csv')
ref_genome_fasta = '/home/toomas/git/RNA-tools/Toomas/genbank_mg1655.fasta'

selected_column = 'rrsA'
selected_positions_column = gene_dic[selected_column].column_label

def by_strand_sorter(index, add_to_pos_reads, add_to_neg_reads, counter_pos, counter_neg):
    if index in counter_pos and index in counter_neg:
        add_to_pos_reads.append(counter_pos[index])
        add_to_neg_reads.append(counter_neg[index])
    elif index not in counter_pos and index in counter_neg:
        add_to_neg_reads.append(counter_neg[index])
        add_to_pos_reads.append(1)
    elif index in counter_pos and index not in counter_neg:
        add_to_pos_reads.append(counter_pos[index])
        add_to_neg_reads.append(1)
    else:
        add_to_pos_reads.append(1)
        add_to_neg_reads.append(1)

def data_generator(index, gene_dic, selected_column, selected_positions_column):
    
    #Creates dataframe from csv
    df = pd.read_csv(index)
    
    #Fill the blank spaces with 0 in dataframe
    selected_df = df.fillna(0)
    
    #Create a new column in dataframe, which contains the number of rRNA operons which read maps to. 
    df['gene_count_str'] = df['rrsA'] + df['rrsB'] + df['rrsC'] + df['rrsD'] + df['rrsE'] + df['rrsG'] + df['rrsH']
    #Creates a new fataframe which contains only reads which map to only a certain numbero f rRNA operons
    df = df.query('gene_count_str == 7')
    
    values_list_pos = []
    values_list_neg = []
    for value in selected_df[selected_positions_column].dropna():
        if '+AC0-' in str(value):                
            value = str(value)[4:]
            values_list_pos.append(int(float(value)))
        elif '-' in str(value):
            value = str(value)[1:]
            values_list_neg.append(int(float(value)))
        else:
            values_list_pos.append(int(float(value)))
    
    counter_pos = Counter(values_list_pos)
    counter_neg = Counter(values_list_neg)

    #Create empty lists for storing reads 5 prime position count on genome.
    # nucl_data = genome position, y = count, y_read_count = total count for both cutting site and other positions
    #cut = reads starting with cutting site, pos = positive strand, neg = negative strand
    nucl_data = []
    y_pos_read_count = []
    y_pos = []
    y_pos_cut = []
    y_neg_read_count = []
    y_neg = []
    y_neg_cut = []
    
    for fasta in SeqIO.parse(ref_genome_fasta, "fasta"):
        pass
    
    for index in range(gene_dic[selected_column].start_pos, gene_dic[selected_column].end_pos):

        if str(fasta.seq[index-1: index+2]) == 'ACA':
            y_pos.append(1)
            y_neg.append(1)
            by_strand_sorter(index, y_pos_cut, y_neg_cut, counter_pos, counter_neg) 
        else:
            y_pos_cut.append(1)
            y_neg_cut.append(1)
            by_strand_sorter(index, y_pos, y_neg, counter_pos, counter_neg)
            
        nucl_data.append(index)
        by_strand_sorter(index, y_pos_read_count, y_neg_read_count, counter_pos, counter_neg)
    
    return (nucl_data, y_pos_read_count ,y_pos, y_pos_cut, y_neg_read_count, y_neg, y_neg_cut)



nucl_data_1, y_pos_read_count_1, y_pos_1, y_pos_cut_1, y_neg_read_count_1, y_neg_1, y_neg_cut_1 = data_generator(df1, gene_dic, 
                                                                        selected_column, 
                                                                        selected_positions_column)

nucl_data_2, y_pos_read_count_2, y_pos_2, y_pos_cut_2, y_neg_read_count_2, y_neg_2, y_neg_cut_2 = data_generator(df2, gene_dic, 
                                                                     selected_column, 
                                                                     selected_positions_column)

MA_X = [(math.log(float(y1), 2) + math.log(float(y2), 2))/2 for y1, y2 in zip(y_pos_1, y_pos_2)]
MA_Y = [math.log((float(y1)/y2), 2) for y1, y2 in zip(y_pos_1, y_pos_2)]
MA_cut_X = [(math.log(float(y1), 2) + math.log(float(y2), 2))/2 for y1, y2 in zip(y_pos_cut_1, y_pos_cut_2)]
MA_cut_Y= [math.log((float(y1)/y2), 2) for y1, y2 in zip(y_pos_cut_1, y_pos_cut_2)]

def onpick3(event):
    ind = event.ind
    nucleotide_pos = np.take(nucl_data_1, ind)
    sample_1_read_count = np.take(y_pos_read_count_1, ind)
    sample_2_read_count = np.take(y_pos_read_count_2, ind)
    for array_ind in range(len(nucleotide_pos)):
        print "Nucleotide position: {0} , s1 count: {1} , s2 count {2}"\
                .format(nucleotide_pos[array_ind], sample_1_read_count[array_ind], 
                sample_2_read_count[array_ind])

fig, ax = plt.subplots()
ax.set_ylabel('M', fontsize=25)
ax.set_xlabel('A', fontsize=25)
ax.scatter(MA_X, MA_Y, alpha=0.5, linewidths=( 0, 0, 0), picker=True)
ax.scatter(MA_cut_X, MA_cut_Y, c='r', alpha=0.5, linewidths=( 0, 0, 0), picker=True)
fig.canvas.mpl_connect('pick_event', onpick3)

fig.savefig('test5.pdf')
plt.show()