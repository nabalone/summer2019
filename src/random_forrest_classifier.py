# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 12:51:52 2019

@author: Noel
"""

import csv
import os
from sklearn import svm
from sklearn.metrics import confusion_matrix
from sklearn.decomposition import PCA
import numpy as np
import matplotlib.pyplot as plt
from sklearn.utils.multiclass import unique_labels
from sklearn.model_selection import LeaveOneOut, cross_val_predict, StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
import random
from imblearn.over_sampling import SMOTE
from sklearn.model_selection import LeaveOneOut,train_test_split
from sklearn.ensemble import RandomForestClassifier
import scipy
import math
from sklearn import preprocessing
import random
import pandas as pd

#from plot_cm import *

CLASS_WEIGHT = 'balanced'
WHITEN = True
USE_RF = False #otherwise SVM
PROJ_HOME = os.environ['DATA_SRCDIR']
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
 'KronMag_3', 'Abs. Mag_3','Ellipticity_3', 'Z_3', 'pixelRank_3', 'chanceCoincidence_3']
cols = cols1[:]
for i in range(4,7):
    for e in cols1:
        cols.append(e[:-1] + str(i))

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
    print(data)
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

'''made up metric for loosely ranking confusion matrices'''
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
    print(X.shape)
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
        loo = LeaveOneOut()
        skf = StratifiedKFold(n_splits=9)

        y_pred = np.zeros(len(y)) - 1
        count=0
        for train_index, test_index in loo.split(X, y):
            print(count)
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
            
        np.save(name_extension + '_ytrue', y)
        np.save(name_extension + '_ypred', y_pred)
        #plot_confusion_matrix(y, y_pred, name_extension=name_extension)
        importances = clf.feature_importances_
        print(importances)
        print(y_pred)
        print(y_test)
        plt.bar(range(len(importances)), importances)
        print(importances)
        importances = np.array(importances)
        np.save(DESTDIR + "/importances" + ext, importances, allow_pickle=True, fix_imports=True)
        plt.savefig(DESTDIR + "/importances_final%s.png" % ext)
        print(y)
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
    
    run(X1, y_ia, 200, DESTDIR + '/ia' + ext)

    end = time.time()
    print(end - start)

if __name__=='__main__':
    main()
