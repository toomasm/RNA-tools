import matplotlib
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt

plt.plot([1,2,3,4], [1,4,9,16], 'r-',  linewidth=2.0)
plt.axis([0, 6, 0, 20])
plt.show()
