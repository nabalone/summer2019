import matplotlib.pyplot as plt
import os
import glob
import numpy as np

def make_plot(filename, stack=False):
    print(filename)
    if stack:
        fils = glob.glob(filename[:-5] + '*' + filename[-4:])
        if len(fils) != 4:
            raise
        vals = []
        accs = []
        for my_filename in fils:
            my_vals = []
            my_accs = []
            with open(my_filename, 'r') as f:
                for line in f:
                    if 'val_acc' in line:
                        acc_ind = line.index('- acc')
                        acc = float(line[acc_ind+7:acc_ind+13])
                        my_accs.append(acc)
                        val_ind = line.index('val_acc')
                        val = float(line[val_ind+9:val_ind+15])
                        my_vals.append(val)
            print(len(vals))
            print(len(my_vals))
            if len(vals) != len(my_vals):
                print(my_vals[:3])
                print(my_vals[-3:])
            if len(vals) == 0:
                vals = np.array(my_vals)
                accs = np.array(my_accs)
            else:
                accs = accs + np.array(my_accs)
                vals = vals + np.array(my_vals)
                
    else:            
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

    plt.plot(1-np.array(accs)/5, 'ro')
    plt.plot(1-np.array(vals)/5, 'bo')
    plt.axis([0, 300, 0, 1])
    plt.grid(which='both')
    plt.savefig(filename[:-4] + '.png')#[-9:-4]+'_3')
    print(filename[:-4] + '.png')
    plt.close()
    
#for letter in ['a','b','c','d','e','f','g','h','i','j']:
    #make_plot(os.getcwd() + '/second_run/cnn_%s.log' % letter)
    
files = glob.glob('sixth_run/cnn_run*0.log')
count = 0
for fil in files:
    make_plot(fil, stack=True)
    count += 1
    print(count)


