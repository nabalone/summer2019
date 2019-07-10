# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 12:51:52 2019

@author: Faith
"""
import csv
from sklearn import svm
from sklearn.metrics import confusion_matrix
from sklearn.decomposition import PCA

headers = 'ID,KronRad_3,separation_3,x_3,y_3,KronMag_3,Angle_3,Ellipticity_3,RA_3,Host RA_3,DEC_3,Host Dec_3,Discrepency_3,Z_3,Hecto Z_3,KronRad_4,separation_4,x_4,y_4,KronMag_4,Angle_4,Ellipticity_4,RA_4,Host RA_4,DEC_4,Host Dec_4,Discrepency_4,Z_4,Hecto Z_4,KronRad_5,separation_5,x_5,y_5,KronMag_5,Angle_5,Ellipticity_5,RA_5,Host RA_5,DEC_5,Host Dec_5,Discrepency_5,Z_5,Hecto Z_5,KronRad_6,separation_6,x_6,y_6,KronMag_6,Angle_6,Ellipticity_6,RA_6,Host RA_6,DEC_6,Host Dec_6,Discrepency_6,Z_6,Hecto Z_6'
headers = headers.split(',')

def chooseX():
    

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


    
import numpy as np
import matplotlib.pyplot as plt

from sklearn.utils.multiclass import unique_labels

type_to_int = {'SNIa':0, 'SNIbc':1,'SNII':2, 'SNIIn':3,  'SLSNe':4}
    
def plot_confusion_matrix(y_true, y_pred, cmap=plt.cm.Blues):
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
    print(cm)
    
    fig, ax = plt.subplots()
    im = ax.imshow(cm, interpolation='nearest', cmap=cmap)
    ax.figure.colorbar(im, ax=ax)
    # We want to show all ticks...
    ax.set(xticks=np.arange(cm.shape[1]),
           yticks=np.arange(cm.shape[0]),
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
    return ax
    

def pca_whiten(X):
    X = np.array(X)
    return PCA(whiten=True).fit_transform(X)        
            



X = []
y = []
X2 = []
y2 = []
with open('goodSecondRun/galaxiesdata.csv') as csvfile:
    csvreader = csv.reader(csvfile)
    #throw away header and blank line
    csvreader.next()
    csvreader.next()
    for i in range(250):
        row = csvreader.next()
        for j in range(len(row)):
            try:
                row[j] = float(row[j])
            except:
                row[j] = 0.
        if len(row) == 57:
            X.append([row[13]])
            y.append(type_to_int[typeDict[pad(int(row[0]))]])
        
    clf = svm.SVC(gamma='scale', class_weight = 'balanced')#{4: 1., 2: 1., 3: 1., 0: 0.222, 1: 2.5})
#{4: 1/13., 2: 1/93., 3: 1/25., 0: 1/357., 1: 1/18.}
    X = pca_whiten(X)
    clf.fit(X, y) 
    
    for i in range(249):
        row = csvreader.next()
        for j in range(len(row)):
            try:
                row[j] = float(row[j])
            except:
                row[j] = 0.
        if len(row) == 57:
            X2.append([row[13]])
            y2.append(type_to_int[typeDict[pad(int(row[0]))]])
    X2 = pca_whiten(X2)        
    print clf.score(X2, y2)
    
    y2p = clf.predict(X2)
    print plot_confusion_matrix(y2, y2p)
    
   
    

