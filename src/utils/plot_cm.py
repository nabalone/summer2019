import os
import glob
from sklearn.metrics import confusion_matrix
import numpy as np
import matplotlib.pyplot as plt

PROJ_HOME = os.environ['DATA_SRCDIR']

SOURCE_DIR = PROJ_HOME + '/src/outputs/cnn_kfold_results/'
SOURCE_DIR_IA = PROJ_HOME + '/src/outputs/cnn_kfold_results_ia/'
OUTPUT_DIR = PROJ_HOME + '/src/outputs/'

ADD_RANDOM = 0

def pad(n):
    n = str(n)
    while len(n) < 6:
        n = '0' + n
    return n

def diagonalishness(m):
    count = 0
    for i in range(len(m)):
        if np.argmax(m[:,i]) == i:
            count += 1
        if np.argmax(m[i,:]) == i:
            count += 1
    return count

def balanced_score(m):
    count = 0
    for i in range(len(m)):
        count += m[i][i]
    return count/len(m)

namecount = 0
def namegen():
    global namecount
    namecount += 1
    return str(namecount)
        
def parser(filename, num, stack=False):
    if stack:
        cm = []#np.zeros(5,5)
        fils = glob.glob(filename[:-5] + '*' + filename[-4:])
        if len(fils) != 4:
            #TODO RESTORE
            pass#raise Exception
        for my_filename in fils:
            fil = open(my_filename, 'r')
            print(num)
            print(type(num))
            lines = fil.readlines()[-1*num:]
            my_cm = []
            for line in lines:
                row = []
                cut_line = line[2:-2]
                cut_line = cut_line.split(']')[0]
                print(cut_line)
                nums = cut_line.split()
                for i in nums:
                    row.append(int(i))
                    
                my_cm.append(row)
        #cms.append(np.array(my_cm))
            if len(cm)==0:
                cm = np.array(my_cm)
            else:
                    cm = cm + np.array(my_cm)
            
    else:       
        fil = open(filename, 'r')
        lines = fil.readlines()[-1*num:]
        cm = []
        for line in lines:
            row = []
            cut_line = line[2:-2]
            cut_line = cut_line.split(']')[0]
            print(cut_line)
            nums = cut_line.split()
            for i in nums:
                row.append(int(i))
            cm.append(row)
            
    return(np.array(cm))
    

def plot_confusion_matrix(y_true, y_pred, cmx = None, to_include='', name_extension='', 
                          cmap=plt.cm.Blues, importances=None, parsed=None):
    """
    This function prints and plots the confusion matrix.
    Normalization can be applied by setting `normalize=True`.
    """
    if str(parsed) != 'None':
        cm = parsed

    else:
        # Compute confusion matrix
        cm = confusion_matrix(y_true, y_pred)# if not cmx.any else cmx
    cm = np.array(cm)
    #normalize
    cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
    if len(cm) == 5:
        classes = np.array(['SNIa', 'SNIbc','SNII', 'SNIIn',  'SLSNe'])
        order = [4,2,3,0,1]
    elif len(cm) == 4:
        classes = np.array(['SNIa', 'SNIbc','SNII', 'SNIIn'])
        order = [2,3,0,1]
    elif len(cm) == 3:
        classes = np.array(['SNIa', 'Ibc/SLS','II/IIn'])
        order = [1,2,0]
    elif len(cm) ==2:
        order = [1,0]
        classes = np.array(['SNIa','CC'])
    else:
        raise

    classes = classes[order]
    print('Confusion matrix:')
    print(cm)

    bal_score = str(balanced_score(cm))[:4]
    diag = str(diagonalishness(cm))
    info_str = "bal_score:" + bal_score \
                + " diag:" + diag + '\n' \
                + str(to_include)
    
    fig, ax = plt.subplots()
    im = ax.imshow(cm, interpolation='nearest', cmap=cmap, vmin=0, vmax=1)

    ax.figure.colorbar(im, ax=ax)
    # We want to show all ticks...
    ax.set(xticks=np.arange(cm.shape[1]),
           yticks=np.arange(cm.shape[0]),
           title=info_str,
           xticklabels=classes, yticklabels=classes,
           ylabel='True label',
           xlabel='Predicted label')
    ax.set_ylim(len(cm) - 0.5, -0.5)

    #ax.tick_params(axis='both', which='major', labelsize=10)
    
    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")
    font = {
    'weight' : 'normal',
    'size'   : 20}
    plt.rc('font', **font)
    # Loop over data dimensions and create text annotations.
    fmt = '.2f' if True else 'd'
    thresh = cm.max() / 2.
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            ax.text(j, i, format(cm[i, j], fmt),
                    ha="center", va="center",
                    color="white" if cm[i, j] > thresh else "black")
    
    #fig.tight_layout()
    
    plt.savefig(name_extension +'.png', bbox_inches = "tight")
    plt.grid(False)
    plt.show()
    plt.close()    
    
def main():
    def combine(files):
        all_pred = []
        all_true = []
        for predfil in files:
            truefil = predfil.replace('pred','true')
            pred = list(np.argmax(np.load(predfil), axis=1))
            true = list(np.load(truefil))
            all_pred.extend(pred)
            all_true.extend(true)
        return(all_true, all_pred)
            
    files_ia = glob.glob(SOURCE_DIR_IA + 'y_pred_ia_fold*.npy')
    files = glob.glob(SOURCE_DIR + 'y_pred_fold*.npy')
    if files_ia:
        print(1)
        true_ia, pred_ia = combine(files_ia)
        plot_confusion_matrix(true_ia, pred_ia, 
                              name_extension=(OUTPUT_DIR + 'cnn_cm_ia'))
    if files:
        print(2)
        true_5, pred_5 = combine(files)
        plot_confusion_matrix(true_5, pred_5, name_extension=OUTPUT_DIR + 'cnn_cm')

    
if __name__ == '__main__':
    main()