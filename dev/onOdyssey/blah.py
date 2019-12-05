import matplotlib.pyplot as plt
import os
import glob
import numpy as np

def make_plot(filename):
    vals = []
    accs = []
    with open(filename, 'r') as f:
        for line in f:
            if 'val_acc' in line:
                acc_ind = line.index('- acc')
                acc = float(line[acc_ind+7:acc_ind+13])
                accs.append(acc)
                val_ind = line.index('val_acc')
                val = float(line[val_ind+9:val_ind+15])
                vals.append(val)
    plt.plot(1-np.array(accs), 'ro')
    plt.plot(1-np.array(vals), 'bo')
    plt.axis([0, 200, 0, 1])
    plt.grid(which='both')
    plt.savefig(filename[:-4] + '.png')#[-9:-4]+'_3')
    print(filename[:-4] + '.png')
    plt.close()
    
#for letter in ['a','b','c','d','e','f','g','h','i','j']:
    #make_plot(os.getcwd() + '/second_run/cnn_%s.log' % letter)
    
files = glob.glob('fifth_run/cnn_run*.log')
count = 0
for fil in files:
    make_plot(fil)
    count += 1
    print(count)


