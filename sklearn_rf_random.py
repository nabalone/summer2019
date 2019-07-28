# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 12:51:52 2019

@author: Faith
"""
#to_include = ['Z', 'Abs. Mag', 'separation (kpc)', 'area (kpc^2)', 'KronRad (kpc)']
#fig_name = "redshiftKronmagKronradSeparation_svm_weightmess_4th"
CLASS_WEIGHT = 'balanced' #{0:0.0000045,  1:5540000000,  2:111,  3:1000000000,  4:768}) FOR RF
        #{0:0.29,  1:5.54,  2:1.11,  3:4.34,  4:7.68}
        #balanced: n_samples / (n_classes * np.bincount(y))
        #0:0.28,  1:5.44,  2:1.07,  3:3.92,  4:7.54
WHITEN = True
USE_RF = False #otherwise SVM
CSV_FILE = '/goodFifthRun/galaxiesdata.csv'
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
from imblearn.over_sampling import SMOTE
from sklearn.model_selection import LeaveOneOut,train_test_split
from sklearn.ensemble import RandomForestClassifier
import scipy
import math
from sklearn import preprocessing
import random 
import pandas as pd

random.seed(7)
PLOT_DIR = os.getcwd() + '/confusions_RF_whitened/'
if not os.path.isdir(PLOT_DIR):
    os.mkdir(PLOT_DIR)

'''unused'''
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
cols1 = ['KronRad (kpc)_3', 'separation (kpc)_3', 'area (kpc^2)_3', 'sep/sqrt(area) (kpc)_3', \
 'KronMag_3', 'Abs. Mag_3','Ellipticity_3', 'Z_3', 'pixelRank_3', 'chance coincidence_3']
cols = cols1[:]
for i in range(4,7):
    for e in cols1:
        cols.append(e[:-1] + str(i))
        
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
    
'''make event id number in int form into string form, 0 padded'''
def pad(n):
    n = str(n)
    while len(n) < 6:
        n = '0' + n
    return n

'''unused'''
def chooseProps(row, to_include):
    limited_row = []
    for prop in to_include:
        if type(prop) == int:
            prop = headers[prop]
        for ind in indexDict[prop]:
            limited_row.append(row[ind])
    return limited_row

'''unused'''
def chooseAll_bad(row, num_random):
    limited_row = []
    x = 19
    for i in range(4):
        offset = 1 + x*i
        limited_row.extend(row[offset: offset + 4])
        limited_row.extend(row[offset + 6: offset + 8])
        limited_row.append(row[offset + 9])
        #limited_row.append(row[offset + 15])
        limited_row.append(row[offset + 17])
        limited_row.append(row[offset + 18])
    for i in range(num_random):
        limited_row.append(random.random())
    return limited_row

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
    #add num_random random numbers to end
    for i in range(num_random)
        rand = np.random.normal(size=len(y))
        X = np.vstack((X.T, rand)).T
    return (X, y)
    

type_to_int = {'SNIa':0, 'SNIbc':1,'SNII':2, 'SNIIn':3,  'SLSNe':4}

'''made up metric for loosely ranking confusion matrices'''
def diagonalishness(m):
    count = 0
    for i in range(len(m)):
        if np.argmax(m[:,i]) == i:
            count += 1
        if np.argmax(m[i,:]) == i:
            count += 1
    return count

'''correctness per type average'''
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


def plot_confusion_matrix(y_true, y_pred, name_extension, cmap=plt.cm.Blues):
    """
    This function prints and plots the confusion matrix.
    Normalization can be applied by setting `normalize=True`.
    """

    # Compute confusion matrix
    cm = confusion_matrix(y_true, y_pred)

    #normalize
    cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
    print(str(cm))
    classes = np.array(['SNIa', 'SNIbc','SNII', 'SNIIn',  'SLSNe'])
    classes = classes[unique_labels(y_true, y_pred)]

    bal_score = str(balanced_score(cm))[:4]
    diag = str(diagonalishness(cm))
    info_str = "bal_score:" + bal_score \
                + " diag:" + diag + '\n' \
                + str(name_extension)
    print(info_str)
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
    plt.savefig(PLOT_DIR + name_extension + '.png')
    plt.close()
#    to_write = [number, bal_score, diag]
#    included_set = set(to_include)
#    for h in headers:
#        if h in included_set:
#            to_write.append('x')
#        else:
#            to_write.append(None)
#    csvwriter.writerow(to_write)
#    return ax


# def resample(X, y):
#     sm = SMOTE(sampling_strategy='not majority', random_state=7)
#     X_res, y_res = sm.fit_resample(X, y)
#     random.seed(7)
#     random.shuffle(X_res)
#     random.seed(7)
#     random.shuffle(y_res)
#     return (X_res, y_res)

'''unused'''
def pca_whiten(X):
    X = np.array(X)
    return PCA(whiten=True).fit_transform(X)

def collect(num_random):
    X = []
    y = []

    '''
    with open(os.getcwd() + '/goodSeventhRun/galaxiesdata7_no_outliers.csv') as csvfile:
        csvreader = csv.reader(csvfile)
        #throw away header and blank line
        csvreader.__next__()
        #print csvreader.next()
        for row in csvreader:
            # convert all to floats
            for j in range(len(row)):
                try:
                    row[j] = float(row[j])
                    if np.isnan(row[j]):
                        row[j] = 0.
                except:
                    row[j] = 0
            # for now only use complete rows
            if len(row) == 1 + 4*len(headers):
                X.append(chooseAll(row, num_random))
                y.append(type_to_int[typeDict[pad(int(row[0]))]])
    #X = pca_whiten(X)
    '''
    X, y = chooseAll(os.getcwd() + CSV_FILE, num_random)
    print(X.shape)
    X = preprocessing.scale(X)
    X = np.array(X)
    y = np.array(y)
    
    #np.random.shuffle(y)
    return (X, y)

def run(X, y, n_est, name_extension):
#       clf = RandomForestClassifier(n_estimators = n_estimators, class_weight = 'balanced_subsample')
        loo = LeaveOneOut()
        skf = StratifiedKFold(n_splits=9)

        y_pred = np.zeros(len(y)) - 1
        count=0
        for train_index, test_index in skf.split(X, y):
            print(count)
            count += 1
            
            #print('Currently training',test_index[0],' of ',len(X))
            X_train, X_test = X[train_index], X[test_index]
            y_train, y_test = y[train_index], y[test_index]

            sampler = SMOTE(sampling_strategy='all')#'not majority') #{0:290,1:10000,2:1000,3:10000,4:1000}, random_state=7)
            #X_res, y_res = X_train, y_train
            X_res, y_res = sampler.fit_resample(X_train, y_train)
            new_ind = np.random.permutation(range(len(y_res)))
            X_res = X_res[new_ind,:]
            y_res = y_res[new_ind]
            if count < 2:
                np.save("deletable_yres", np.array(y_res))
            clf = RandomForestClassifier(n_estimators=n_est)
                    #class_weight = {0:0.00001, 1:2000000000, 2:100, 3:3000000000, 4:2000})
            clf.fit(X_res,y_res)
            y_pred[test_index] = clf.predict(X_test)
            #print(y_pred[test_index])
            #print(y_test)
            #print('\n')
        np.save(name_extension, y_pred)
        plot_confusion_matrix(y, y_pred, name_extension)
        importances = clf.feature_importances_
        print(y_pred)
        plt.bar(range(len(importances)), importances)
        print(importances)
        plt.savefig("importances.png")
        print(y)

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
#print(y1)
# print('c')
# run(X0, y0, 100, 'rf_0rands_100ests_a')
# print('d')
# run(X0, y0, 100, 'rf_0rands_100ests_b')
# print('e')
run(X1, y1, 100, 'rf_1rands_50ests_a')
# print('f')
# run(X1, y1, 100, 'rf_1rands_100ests_b')
# print('a')
# run(X0, y0, 200, 'rf_0rands_200ests_a')
# print('b')
# run(X0, y0, 200, 'rf_0rands_200ests_b')
end = time.time()
print(end - start)
