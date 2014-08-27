import pandas as pd
import pysam
from sets import Set
from progressbar import AnimatedMarker, Bar, BouncingBar, Counter, ETA, \
    ProgressBar, ReverseBar, RotatingMarker, \
    SimpleProgress, Timer


widgets = [Bar('>'), ' ', ETA(), ' ', ReverseBar('<')]
pbar = ProgressBar(widgets=widgets, maxval=2*22464520).start()
i = 0  

samfile = pysam.Samfile('/media/markus/DATA2/MazF_PNK_sorted.bam', 'rb')
names = set()

columns = ['rrsA', 'rrsB', 'rrsC', 'rrsD', 'rrsE', 'rrsG', 'rrsH']

for read in samfile:
    names.add(read.qname)
    i += 1
    pbar.update(i)

     
names_list = list(names)

df = pd.DataFrame(index=names_list, columns=columns)

samfile = pysam.Samfile('/media/markus/DATA2/MazF_PNK_sorted.bam', 'rb')
for read in samfile:
    i += 1
    pbar.update(i)
    if (read.pos > 4035153 and read.pos < 4040906):
        df['rrsA'][read.qname] = 1
    elif (read.pos > 4165428 and read.pos < 4172057):
        df['rrsB'][read.qname] = 1
    elif (read.pos > 3941327 and read.pos < 3946872):
        df['rrsC'][read.qname] = 1
    elif (read.pos > 3423194 and read.pos < 3429236):
        df['rrsD'][read.qname] = 1
    elif (read.pos > 4207532 and read.pos < 4213234):
        df['rrsE'][read.qname] = 1
    elif (read.pos > 2725746 and read.pos < 2731600):
        df['rrsG'][read.qname] = 1
    elif (read.pos > 223408 and read.pos < 229167):
        df['rrsH'][read.qname] = 1
    else:
        continue
pbar.finish()        

samfile.close()
df.to_csv('index.csv',index=True,header=True)        