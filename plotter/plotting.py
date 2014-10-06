import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from collections import namedtuple, OrderedDict, Counter



GeneInfo = namedtuple('GeneInfo', ['start_pos', 'end_pos', 'column_label'])

gene_dic = OrderedDict([('rrsA', GeneInfo(4035153, 4040906, 'P1')), 
                        ('rrsB', GeneInfo(4165428, 4172057, 'P2')),
                        ('rrsC', GeneInfo(3941327, 3946872, 'P3')),
                        ('rrsD', GeneInfo(3423194, 3429236, 'P4')),
                        ('rrsE', GeneInfo(4207532, 4213234, 'P5')),
                        ('rrsG', GeneInfo(2725746, 2731600, 'P6')),
                        ('rrsH', GeneInfo(223408, 229167, 'P7'))])

df = pd.read_csv('/media/toomas/DATA1/aligment_bowtie/MazF_2h/aligned/ok_MazF3prime_PNK-_index.csv')

selected_column = 'rrsA'
selected_positions_column = gene_dic[selected_column].column_label

values_list = []
for value in df[selected_positions_column].dropna():
    if '+AC0-' in str(value):
        value = str(value)[4:]
    elif '-' in str(value):
        value = str(value)[1:]
    values_list.append(int(float(value)))

counter = Counter(values_list)
x = []
y = []
#print counter
ptcnt = 0
for index in range(gene_dic[selected_column].start_pos, gene_dic[selected_column].end_pos):
    if index in counter:
        ptcnt += 1
        #print ptcnt
        y.append(counter[index])
    else:
        y.append(0)
    x.append(index)
#print all(v == 0 for v in y)  

fig = plt.figure(figsize=(8,6))
plt.scatter(x, y)
plt.savefig('test2.pdf')
plt.show()
