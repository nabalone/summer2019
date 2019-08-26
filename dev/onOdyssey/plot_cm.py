import csv
import os
import random
from sklearn import svm
from sklearn.metrics import confusion_matrix
from sklearn.decomposition import PCA
import numpy as np
import matplotlib.pyplot as plt
from sklearn.utils.multiclass import unique_labels
from sklearn.model_selection import LeaveOneOut, cross_val_predict
from sklearn.ensemble import RandomForestClassifier
PLOT_DIR = os.getcwd() + '/cm_plot'#\\deletable\\'#'/confusions_svm_whitened_cv100/'
ADD_RANDOM = 0
if not os.path.isdir(PLOT_DIR):
    os.mkdir(PLOT_DIR)

def pad(n):
    n = str(n)
    while len(n) < 6:
        n = '0' + n
    return n

def chooseProps(row, to_include):
    limited_row = []
    for prop in to_include:
        if type(prop) == int:
            prop = headers[prop]
        for ind in indexDict[prop]:
            limited_row.append(row[ind])
    for i in range(ADD_RANDOM):
        limited_row.append(random.random())
    return limited_row

type_to_int = {'SNIa':0, 'SNIbc':1,'SNII':2, 'SNIIn':3,  'SLSNe':4}

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
        

def plot_confusion_matrix(y_true, y_pred, cmx = None, to_include='', name_extension='', 
                          cmap=plt.cm.Blues, importances=None):
    """
    This function prints and plots the confusion matrix.
    Normalization can be applied by setting `normalize=True`.
    """
    
    # Compute confusion matrix
    cm = confusion_matrix(y_true, y_pred)#s if not cmx.any else cmx
    
    #normalize
    #cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
    #classes = np.array(['SNIa', 'SNIbc','SNII', 'SNIIn',  'SLSNe'])
    #classes = classes[unique_labels(y_true, y_pred)]
    classes = np.array(['SNIa','CC'])
    bal_score = str(balanced_score(cm))[:4]
    diag = str(diagonalishness(cm))
    info_str = "bal_score:" + bal_score \
                + " diag:" + diag + '\n' \
                + str(to_include)
    
    fig, ax = plt.subplots()
    im = ax.imshow(cm, interpolation='nearest', cmap=cmap)
    ax.figure.colorbar(im, ax=ax)
    # We want to show all ticks...
    ax.set(xticks=np.arange(cm.shape[1]),
           yticks=np.arange(cm.shape[0]),
           title=info_str,
           xticklabels=classes, yticklabels=classes,
           ylabel='True label',
           xlabel='Predicted label')
    #ax.tick_params(axis='both', which='major', labelsize=10)
    
    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")
    font = {
    'weight' : 'normal',
    'size'   : 10}
    plt.rc('font', **font)
    # Loop over data dimensions and create text annotations.
    fmt = '.2f' if True else 'd'
    thresh = cm.max() / 2.
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            ax.text(j, i, format(cm[i, j], fmt),
                    ha="center", va="center",
                    color="white" if cm[i, j] > thresh else "black")
    fig.tight_layout()
    number = namegen()
    
    plt.savefig(PLOT_DIR + bal_score + '_' + diag + '_' + str(ADD_RANDOM) \
                + name_extension + '.png')#number + '.png')
    plt.show()
    plt.close()    
    print(1)
  
