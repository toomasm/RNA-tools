import pandas as pd
from collections import namedtuple, OrderedDict, Counter, defaultdict
import plotarray
import os


def get_plot_data(selected_column, index_file):

    GeneInfo = namedtuple('GeneInfo', ['start_pos', 'end_pos', 'column_label'])

    gene_dic = OrderedDict([('rrsA', GeneInfo(4035153, 4040906, 'P1')), 
                            ('rrsB', GeneInfo(4165428, 4172057, 'P2')),
                            ('rrsC', GeneInfo(3941327, 3946872, 'P3')),
                            ('rrsD', GeneInfo(3423194, 3429236, 'P4')),
                            ('rrsE', GeneInfo(4207532, 4213234, 'P5')),
                            ('rrsG', GeneInfo(2725746, 2731600, 'P6')),
                            ('rrsH', GeneInfo(223408, 229167, 'P7'))])


    selected_positions_column = gene_dic[selected_column].column_label

    values_list_pos = []
    values_list_neg = []
    for value in index_file[selected_positions_column].dropna():
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
    plot_data_gen['pos']['x'] = x_pos
    plot_data_gen['pos']['y'] = y_pos
    plot_data_gen['neg']['x'] = x_neg
    plot_data_gen['neg']['y'] = y_neg
    return plot_data_gen

    

if __name__ == '__main__':
    def data_dict():
        return defaultdict(data_dict)

    hdf_data = data_dict()
    indexes = ['Index1.csv', 'Index2.csv']
    operons = ['rrsA', 'rrsB', 'rrsC', 'rrsD', 'rrsE', 'rrsG', 'rrsH']
    for index_path in indexes:
        index_basename = os.path.basename(index_path)
        index_ext = os.path.splitext(index_basename)
        index_name = index_ext[0]
        index = pd.read_csv(index_path)
        for operon in operons:
            hdf_data[index_name][operon] = get_plot_data(operon, index)
    h5_filename = (index_name + '_' + str(len(operons)) + '_.hdf5')
    plotarray.core.save_plot_data(hdf_data, h5_filename)