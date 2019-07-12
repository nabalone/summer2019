# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 12:51:52 2019

@author: Faith
"""
#to_include = ['Z', 'Abs. Mag', 'separation (kpc)', 'area (kpc^2)', 'KronRad (kpc)']
#fig_name = "redshiftKronmagKronradSeparation_svm_weightmess_4th"
CLASS_WEIGHT = "balanced"#{0:0.0000045,  1:5540000000,  2:111,  3:1000000000,  4:768} # FOR RF
        #{0:0.29,  1:5.54,  2:1.11,  3:4.34,  4:7.68}
        #balanced: n_samples / (n_classes * np.bincount(y))
        #0:0.28,  1:5.44,  2:1.07,  3:3.92,  4:7.54
WHITEN = True
USE_RF = True #otherwise SVM

ADD_RANDOM = 0

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
PLOT_DIR = os.getcwd() + '\\deletable\\'#'/confusions_svm_whitened_cv100/'
if not os.path.isdir(PLOT_DIR):
    os.mkdir(PLOT_DIR)

headers = 'KronRad (kpc)_3,separation (kpc)_3,area (kpc^2)_3,sep/area (kpc)_3,x_3,y_3,KronMag_3,Abs. Mag_3,Angle_3,Ellipticity_3,RA_3,Host RA_3,DEC_3,Host Dec_3,Discrepency (arcsecs)_3,Z_3,SDSS Photoz_3,pixelRank_3,chance coincidence_3'
headers = headers.split(',')
for i in range(len(headers)):
    headers[i] = headers[i][:-2]
indexDict = {}
for i in range(len(headers)):
    indexDict[headers[i]] = []
    for j in range(4):
        indexDict[headers[i]].append(1+i+j*len(headers))
        
goodHeaderIndices = [0,1,2,7,9,15,17,18]
allCombos = [[]]
for i in goodHeaderIndices:
    copy = allCombos[:]
    for row in copy:
        allCombos.append(row + [headers[i]])

'''load event type dictionary'''
typeDict = {}
typefile = open('ps1confirmed_only_sne_without_outlier.txt', 'r')
typefile.readline() #get rid of header
for line in typefile:
    parts = line.split() 
    eventId = parts[0][3:]
    eventType = parts[1]
    typeDict[eventId] = eventType
typefile.close()

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
        
    
def plot_confusion_matrix(y_true, y_pred, to_include, name_extension='', 
                          cmap=plt.cm.Blues, importances=None):
    """
    This function prints and plots the confusion matrix.
    Normalization can be applied by setting `normalize=True`.
    """
    
    # Compute confusion matrix
    cm = confusion_matrix(y_true, y_pred)
    
    #normalize
    cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
    classes = np.array(['SNIa', 'SNIbc','SNII', 'SNIIn',  'SLSNe'])
    classes = classes[unique_labels(y_true, y_pred)]
    
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
    
    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")
    
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
    #plt.show()
    plt.savefig(PLOT_DIR + bal_score + '_' + diag + '_' + str(ADD_RANDOM) \
                + name_extension + '.png')#number + '.png')
    plt.close()    
    to_write = [number, bal_score, diag]
    included_set = set(to_include)
    for h in headers:
        if h in included_set:
            to_write.append('x')
        else:
            to_write.append(None)
    to_write.extend(importances)
    csvwriter.writerow(to_write)
    return ax
    
    

def pca_whiten(X):
    X = np.array(X)
    return PCA(whiten=True).fit_transform(X)        
            

def run(to_include, N_ESTIMATORS = None, NAME_EXTENSION = ''):
    X = []
    y = []
    
    with open('goodFourthRun/galaxiesdata.csv') as csvfile:
        csvreader = csv.reader(csvfile)
        #throw away header and blank line
        csvreader.next()
        csvreader.next()
        for row in csvreader: 
            # convert all to floats
            for j in range(len(row)):
                try:
                    row[j] = float(row[j])
                    if np.isnan(row[j]):
                        row[j] = 0.
                except:
                    row[j] = 0.
            # for now only use complete rows
            if len(row) == 1 + 4*len(headers):
                X.append(chooseProps(row, to_include))
                y.append(type_to_int[typeDict[pad(int(row[0]))]])
        if USE_RF:
            clf = RandomForestClassifier(n_estimators=N_ESTIMATORS, class_weight = CLASS_WEIGHT)
        else:
            clf = svm.SVC(gamma='scale', class_weight = CLASS_WEIGHT)#{4: 1., 2: 1., 3: 1., 0: 0.222, 1: 2.5})
    #{4: 1/13., 2: 1/93., 3: 1/25., 0: 1/357., 1: 1/18.}
        if WHITEN: 
            X = pca_whiten(X)
        clf.fit(X, y)     
        loo = LeaveOneOut()
        y_pred = cross_val_predict(clf, X, y, cv=loo)
        I = clf.feature_importances_ if USE_RF else None
        plot_confusion_matrix(y, y_pred, to_include,
                              name_extension = NAME_EXTENSION, 
                              importances = I)
        

import time
start = time.time()
with open(PLOT_DIR + "log.csv", "w+") as destfile:
    csvwriter = csv.writer(destfile)
    x = ['number', 'bal_score', 'diag']
    x.extend(headers)
    csvwriter.writerow(x)
    for combo in [['Z', 'Abs. Mag', 'separation (kpc)', 'area (kpc^2)', \
                   'KronRad (kpc)', 'Ellipticity', 'pixelRank', 'chance coincidence'], ['Z', 'pixelRank']]: #allCombos[1:]:
        run(combo, 10, 'rf10')
        run(combo, 100, '100a')
        run(combo, 100, '100b')
        global ADD_RANDOM
        ADD_RANDOM=1
        run(combo, 100, 'r100a')
        run(combo, 100, 'r100b')
        run(combo, 500, '500a')
        run(combo, 500, '500b')
end = time.time()
print(end - start)
        
    