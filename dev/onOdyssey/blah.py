import matplotlib.pyplot as plt
import os
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
    plt.plot(accs, 'ro')
    plt.plot(vals, 'bo')
    plt.savefig(filename[-9:-4])
    plt.close()
    
for letter in ['a','b','c','d','e','f','g','h','i','j']:
    make_plot(os.getcwd() + '/cnn_%s.log' % letter)


