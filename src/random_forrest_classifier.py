# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 12:51:52 2019

@author: Noel
"""

import os
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from imblearn.over_sampling import SMOTE
from sklearn.model_selection import LeaveOneOut, StratifiedKFold
from sklearn import preprocessing
import random
import pandas as pd
import sys

PROJ_HOME = os.environ['DATA_SRCDIR']
sys.path.append(PROJ_HOME)
from src.utils.plot_cm import pad, plot_confusion_matrix
from src.utils.plot_importances import plot_importances

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-k', '--kfold', type=int, required=False, 
                    help='use kfold cross validation with specified number of \
                    folds instead of leave one out cross validation')
parser.add_argument('--no_redshift', action='store_true', \
                    help='do not include redshift directly in the training data.\
                    it will still be used to calculate absolute quantities')
args = parser.parse_args()
args.kfold=3

NUM_TREES = 700
DESTDIR = PROJ_HOME + "/dev/apple_run"
CSV_FILE = DESTDIR + '/galaxiesdata.csv'

INSIDE_ONLY = False
COLLAPSE = False

ext = '3'
if COLLAPSE:
    ext = ext + '_collapsed'
if INSIDE_ONLY:
    ext = ext + '_inside_only'

random.seed(3)
PLOT_DIR = DESTDIR + '/confusions_RF_whitened/'
if not os.path.isdir(PLOT_DIR):
    os.mkdir(PLOT_DIR)

cols1 = ['KronRad (kpc)_3', 'separation (kpc)_3', 'area (kpc^2)_3', 'sep/sqrt(area) (kpc)_3', \
 'KronMag_3', 'Abs. Mag_3','Ellipticity_3', 'pixelRank_3', 'chanceCoincidence_3']
if args.no_redshift:
    print("Not including redshift.")
else:
    print("Including redshift.")
    cols1.append('redshift_3')

cols = cols1[:]
for i in range(4,7):
    for e in cols1:
        cols.append(e[:-1] + str(i))
       
names = []
for i in cols:
    new=i.replace('3', 'g').replace('4', 'r').replace('5', 'i').replace('6', 'z')
    names.append(new)

'''load event type dictionary'''
typeDict = {}
typefile = open(PROJ_HOME + '/src/ps1confirmed_added.txt', 'r')
typefile.readline() #get rid of header
for line in typefile:
    parts = line.split()
    eventId = parts[0][3:]
    eventType = parts[1]
    typeDict[eventId] = eventType
typefile.close()

if INSIDE_ONLY:
    from inside_and_outside import INSIDE
    ins_set = INSIDE['SNIa'] | INSIDE['SNIbc'] | INSIDE['SNII'] | INSIDE['SNIIn'] 

        


#'''make event id number in int form into string form, 0 padded'''
#def pad(n):
#    n = str(n)
#    while len(n) < 6:
#        n = '0' + n
#    return n

'''getting columns from csv file with pandas'''
def chooseAll(csvfile, num_random):
    data = pd.read_csv(csvfile)
    X = data.loc[:, cols].as_matrix()
    X = np.nan_to_num(X)
    y = []

              
    # lookup identified type of event and add to y
    for i in data['ID']:
        y.append(type_to_int[typeDict[pad(int(i))]])
    np.save(DESTDIR + '/y_actual' + ext, y, allow_pickle=True)
    #add num_random random numbers to end
    for i in range(num_random):
        rand = np.random.normal(size=len(y))
        X = np.vstack((X.T, rand)).T
        
    # if INSIDE ONLY, filter in only events that are inside their host 
    if INSIDE_ONLY:
        new_X = []
        new_y = []
        for i in range(len(X)):
            j = data['ID'][i]
            id_num_str = pad(int(j))
            if id_num_str in ins_set:
                new_X.append(X[i])
                new_y.append(type_to_int[typeDict[id_num_str]])
        X = np.array(new_X)
        y=new_y
        
    return (X, y)


type_to_int = {'SNIa':0, 'SNIbc':1,'SNII':2, 'SNIIn':3,  'SLSNe':4}

'''made-up metric for loosely ranking confusion matrices'''
#def diagonalishness(m):
#    count = 0
#    for i in range(len(m)):
#        if np.argmax(m[:,i]) == i:
#            count += 1
#        if np.argmax(m[i,:]) == i:
#            count += 1
#    return count
#
#'''correctness per type average'''
#def balanced_score(m):
#    count = 0
#    for i in range(len(m)):
#        count += m[i][i]
#    return count/len(m)
#
#namecount = 0
#def namegen():
#    global namecount
#    namecount += 1
#    return str(namecount)


#def plot_confusion_matrix(y_true, y_pred, name_extension, cmap=plt.cm.Blues):
#    # Copied and adapted from
#    #https://scikit-learn.org/stable/auto_examples/model_selection/plot_confusion_matrix.html#sphx-glr-auto-examples-model-selection-plot-confusion-matrix-py
#
#    # Compute confusion matrix
#    cm = confusion_matrix(y_true, y_pred)
#
#    #normalize
#    cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
#    print(str(cm))
#    classes = np.array(['SNIa', 'SNIbc','SNII', 'SNIIn',  'SLSNe'])
#    classes = classes[unique_labels(y_true, y_pred)]
#
#    bal_score = str(balanced_score(cm))[:4]
#    #diag = str(diagonalishness(cm))
#    info_str = "bal_score:" + bal_score
#                #+ " diag:" + diag + '\n'
#                #+ str(name_extension)
#    print(info_str)
#    fig, ax = plt.subplots()
#    im = ax.imshow(cm, interpolation='nearest', cmap=cmap)
#    ax.figure.colorbar(im, ax=ax)
#    # We want to show all ticks...
#    ax.set(xticks=np.arange(cm.shape[1]),
#           yticks=np.arange(cm.shape[0]),
#           title=info_str,
#           xticklabels=classes, yticklabels=classes,
#           ylabel='True label',
#           xlabel='Predicted label')
#
#    # Rotate the tick labels and set their alignment.
#    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
#             rotation_mode="anchor")
#
#    # Loop over data dimensions and create text annotations.
#    fmt = '.2f' if True else 'd'
#    thresh = cm.max() / 2.
#    for i in range(cm.shape[0]):
#        for j in range(cm.shape[1]):
#            ax.text(j, i, format(cm[i, j], fmt),
#                    ha="center", va="center",
#                    color="white" if cm[i, j] > thresh else "black")
#    fig.tight_layout()
#    plt.savefig(PLOT_DIR + 'cm%s.png' % ext)
#    plt.close()

def collect(num_random):
    X = []
    y = []
    X, y = chooseAll(CSV_FILE, num_random)
    #print(X.shape)
    X = preprocessing.scale(X)
    X = np.array(X)

    if COLLAPSE:
        for i in range(len(y)):
            if y[i] > 0:
                y[i]=1

    y = np.array(y)

    return (X, y)


'''Evaulates classifier using Leave One Out cross-validation.
Generates confusion matrix (visualization of success and error) and plot of
feature importances to classification'''
#and importance plots
def run(X, y, n_est, name_extension):
    # with help from Victory Ashley Villar
    
        if args.kfold:
            validator = StratifiedKFold(n_splits=args.kfold)
        else:
            validator = LeaveOneOut()
        
        y_pred = np.zeros(len(y)) - 1
        y_true = np.zeros(len(y)) - 1
        count=0
        for train_index, test_index in validator.split(X, y):
            print("fold: %s"%count)
            count += 1

            X_train, X_test = X[train_index], X[test_index]
            y_train, y_test = y[train_index], y[test_index]

            sampler = SMOTE(sampling_strategy='all')
            X_res, y_res = sampler.fit_resample(X_train, y_train)
            new_ind = np.random.permutation(range(len(y_res)))
            X_res = X_res[new_ind,:]
            y_res = y_res[new_ind]
            if count < 2:
                np.save(DESTDIR + "/yres", np.array(y_res))
            clf = RandomForestClassifier(n_estimators=n_est)
            clf.fit(X_res,y_res)
            y_pred[test_index] = clf.predict(X_test)
            y_true[test_index] = y_test
            
        np.save(name_extension + '_ytrue', y_true)
        np.save(name_extension + '_ypred', y_pred)
        plot_confusion_matrix(y, y_pred, name_extension=name_extension)
        importances = clf.feature_importances_
        #plt.bar(range(len(importances)), importances, color='m')
        importances = np.array(importances)
        plot_importances(importances, names, DESTDIR + "/importances%s.png" % ext)
        np.save(DESTDIR + "/importances" + ext, importances, allow_pickle=True, fix_imports=True)
def main():
    import time
    start = time.time()
    #with open(PLOT_DIR + "log.csv", "w+") as destfile:
    #    csvwriter = csv.writer(destfile)
    #    x = ['number', 'bal_score', 'diag']
    #    x.extend(headers)
    #    csvwriter.writerow(x)
    #    for combo in allCombos[1:]:
    #        run(combo)
    #X0, y0 = collect(0)
    X1, y1 = collect(1)

    #run(X1, y1, 200, DESTDIR + '/blah' + ext)
    
    y1 = np.array(y1)
    y_ia = np.where(y1==0,0,1)
    
    run(X1, y_ia, NUM_TREES, DESTDIR + '/ia' + ext)

    end = time.time()
    print('running time: %s' % (end - start))

if __name__=='__main__':
    main()
