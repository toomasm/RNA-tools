import matplotlib
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
import numpy as np
import pandas as pd
from collections import namedtuple, OrderedDict, Counter, defaultdict
import plotarray


def get_plot_data():

    GeneInfo = namedtuple('GeneInfo', ['start_pos', 'end_pos', 'column_label'])

    gene_dic = OrderedDict([('rrsA', GeneInfo(4035153, 4040906, 'P1')), 
                            ('rrsB', GeneInfo(4165428, 4172057, 'P2')),
                            ('rrsC', GeneInfo(3941327, 3946872, 'P3')),
                            ('rrsD', GeneInfo(3423194, 3429236, 'P4')),
                            ('rrsE', GeneInfo(4207532, 4213234, 'P5')),
                            ('rrsG', GeneInfo(2725746, 2731600, 'P6')),
                            ('rrsH', GeneInfo(223408, 229167, 'P7'))])

    #df = pd.read_csv('generated/MazF2h_PNK_trimmed_index.csv')
    df = pd.read_csv('Index.csv') 

    selected_column = 'rrsA'
    selected_positions_column = gene_dic[selected_column].column_label

    values_list_pos = []
    values_list_neg = []
    for value in df[selected_positions_column].dropna():
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
    x_pos = []
    y_pos = []
    x_neg = []
    y_neg = []

    def data_dict():
      return defaultdict(data_dict)

    plot_data_gen = data_dict()

    for index in range(gene_dic[selected_column].start_pos, gene_dic[selected_column].end_pos):
        if index in counter_pos and index in counter_neg:
            y_pos.append(counter_pos[index])
            y_neg.append(counter_neg[index])
        elif index not in counter_pos and index in counter_neg:
          y_neg.append(counter_neg[index])
          y_pos.append(0)
        elif index in counter_pos and index not in counter_neg:
          y_pos.append(counter_pos[index])
          y_neg.append(0)
        else:
            y_pos.append(0)
            y_neg.append(0)
        x_pos.append(index)
        x_neg.append(index)
    plot_data_gen[selected_column + '_pos']['x'] = x_pos
    plot_data_gen[selected_column + '_pos']['y'] = y_pos
    plot_data_gen[selected_column + '_neg']['x'] = x_neg
    plot_data_gen[selected_column + '_neg']['y'] = y_neg
    h5_filename = selected_column + '.hdf5'
    plotarray.core.save_plot_data(plot_data_gen, h5_filename)




#get_plot_data()

plotting_data = plotarray.core.retrieve_plot_data('rrsA.hdf5')

fig = plt.figure(figsize=(20,10))
ax1 = fig.add_subplot(111)
ax2=plt.twinx()
ax1.set_yscale('linear')
ax2.set_yscale('linear')
plt.gca().invert_yaxis()
#align_yaxis(ax1, plotting_data['rrsA_pos']['y'][...], ax2, 0)

ax1.scatter(plotting_data['rrsA_pos']['x'][...], plotting_data['rrsA_pos']['y'][...],
            s=10, c='b', marker="s", label='first', edgecolors='none', alpha=0.3)
ax2.scatter(plotting_data['rrsA_neg']['x'][...], plotting_data['rrsA_neg']['y'][...],
            s=10, c='r', marker="o", label='second', edgecolors='none', alpha=0.3)
plt.legend(loc='upper left');

axcolor = 'lightgoldenrodyellow'

rax = plt.axes([0.025, 0.5, 0.15, 0.15], axisbg=axcolor)
radio = RadioButtons(rax, ('log', 'linear'), active=0)
def colorfunc(label):
    ax1.set_yscale(label)
    plt.draw()
radio.on_clicked(colorfunc)
#ax1.figure.show()
#plt.savefig('test2.png')

plt.show()