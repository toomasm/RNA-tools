import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import plotarray

plotting_data = plotarray.core.retrieve_plot_data('Index1_mappedto3_.hdf5')


fig = plt.figure(figsize=(20,10))
ax1 = fig.add_subplot(111)
ax2=plt.twinx()
ax1.set_yscale('log')
ax2.set_yscale('linear')
plt.gca().invert_yaxis()

y_pos = plotting_data['Index1']['mappedto3']['pos']['y'][...]# - plotting_data['Index2']['rrsA']['pos']['y'][...]
y_neg = plotting_data['Index1']['mappedto3']['neg']['y'][...]# - plotting_data['Index2']['rrsA']['neg']['y'][...]

ax1.scatter(plotting_data['Index1']['mappedto3']['pos']['x'][...], y_pos ,
            s=10, c='b', marker="s", label='positive', edgecolors='none', alpha=0.3)
ax2.scatter(plotting_data['Index1']['mappedto3']['neg']['x'][...], y_neg,
            s=10, c='r', marker="o", label='negative', edgecolors='none', alpha=0.3)
plt.legend(loc='upper left');
#ax1.figure.show()
plt.savefig('test1.png')

#plt.show()
"""
from pandas import HDFStore

store = HDFStore('Index2_7_.hdf5')

store['table1Name'].to_csv('outputFileForTable1.csv')"""