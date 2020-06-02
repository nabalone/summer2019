# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 12:51:52 2019

@author: Noel
"""

# Note: filters 3,4,5,6 are g,r,i,z respectively
# Types 0,1,2,3,4 are types Ia, Ibc, II, IIn, and superluminous respectively

import os
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from imblearn.over_sampling import SMOTE
from sklearn.model_selection import LeaveOneOut, StratifiedKFold
from sklearn import preprocessing
import random
import pandas as pd
import sys
import time
import pickle

PROJ_HOME = os.environ['DATA_SRCDIR']
sys.path.append(PROJ_HOME)
from src.utils.plot_cm import pad, plot_confusion_matrix
from src.utils.plot_importances import plot_importances
from src.host_property_extractor import COLUMNS

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-k', '--kfold', type=int, required=False, 
                    help='use kfold cross validation with specified number of \
                    folds')
parser.add_argument('--loo', action='store_true', 
                    help='use leave one out cross validation.')
parser.add_argument('--all', action='store_true', 
                    help='train on entire and save model. No validation or testing.')
parser.add_argument('--no_redshift', action='store_true', 
                    help='if this flag is present, we will NOT include redshift \
                    directly in the training data.\
                    it will still be used to calculate absolute quantities')
parser.add_argument('-a', '--ia_only', action='store_true', 
                    help="classify between Type Ia vs. other only")
parser.add_argument('--inside_only', action='store_true', 
                    help="Note this functionality is HACKY and FOR DIAGNOSTIC \
                    PURPOSES ONLY. Only include samples in which SN occured \
                    within host galaxy, to eliminate uncertainty from \
                    calculated the most likely host. Relies on a pre-generated \
                    list of which SN are inside host; I don't remember how I \
                    generated that list")

args = parser.parse_args()


    

NUM_TREES = 700
DESTDIR = PROJ_HOME + "/src/outputs"
CSV_FILE = DESTDIR + '/galaxiesdata.csv'

ext = 'rf'
if args.ia_only:
    ext = ext + '_ia'
if args.inside_only:
    ext = ext + '_inside_only'
if args.kfold:
    ext = ext + '_%skfolds' % args.kfold
else:
    ext = ext + '_loo'

random.seed(3)
PLOT_DIR = DESTDIR + '/confusions_RF_whitened/'
if not os.path.isdir(PLOT_DIR):
    os.mkdir(PLOT_DIR)
    
EXCLUDED_COLS = {'ID', 'comparison', 'RA', 'DEC', 'x', 'y'}
if args.no_redshift:
    EXCLUDED_COLS.add('Redshift')
cols = []
for column in COLUMNS:
    exclude_this_column=False
    for excluded_column in EXCLUDED_COLS:
        if excluded_column in column:
            exclude_this_column=True
            break
    if not exclude_this_column:
        cols.append(column)
print("columns to use: %s" % str(cols))


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

type_to_int = {'SNIa':0, 'SNIbc':1,'SNII':2, 'SNIIn':3,  'SLSNe':4}

if args.inside_only:
    from dev.inside_and_outside import INSIDE
    ins_set = INSIDE['SNIa'] | INSIDE['SNIbc'] | INSIDE['SNII'] | INSIDE['SNIIn'] 

'''getting columns from csv file with pandas'''
def chooseAll(csvfile, num_random, include_ID=False):
    if include_ID:
        cols_to_use = ['ID'] + cols
    else:
        cols_to_use= cols
    data = pd.read_csv(csvfile)
    X_orig = data.loc[:, cols_to_use].values
    X_orig = X_orig.astype(float)
    
    # since previous versions of numpy nan_to_num do not suport posinf args, 
    # cap all values at 1000000. Otherwise we hit errors later
    X_orig = np.where(X_orig>1000000, 1000000, X_orig)
    X_orig = np.where(X_orig<-1000000, -1000000, X_orig)
    #X_orig = data.values
    X_orig = np.nan_to_num(X_orig)#, posinf=1000000., neginf=-1000000.)
    X_orig = np.where(X_orig > 9000000, 9000000, X_orig)
    X_orig = np.where(X_orig > 9000000, 9000000, X_orig)
    y = []
    
    X = X_orig
#TODO this doesn't work because it replaces ID numbers. Find a better solution    
    #there shouldn't be any number greater than about 2000
    #replace any erroneously excessively large numbers with 0
    #X = np.where(X_orig > 100000, 0, X_orig)
#    X = []
#    for i in X_orig:
#        if np.max(i) < 100000:
#            X.append(i)
#    X = np.array(X)
              
    # lookup identified type of event and add to y
    for i in data['ID']:
        y.append(type_to_int[typeDict[pad(int(i))]])
    np.save(DESTDIR + '/y_actual' + ext, y, allow_pickle=True)
    #add num_random random numbers to end
    for i in range(num_random):
        rand = np.random.normal(size=len(y))
        X = np.vstack((X.T, rand)).T
        
    # if INSIDE ONLY, filter in only events that are inside their host 
    if args.inside_only:
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

def collect(num_random):
    X = []
    y = []
    X, y = chooseAll(CSV_FILE, num_random)
    #print(X.shape)
    X = preprocessing.scale(X)
    X = np.array(X)

    if args.ia_only:
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
        if args.kfold:
            val_folds = StratifiedKFold(n_splits=args.kfold).split(X, y)
        elif args.loo:
            val_folds = LeaveOneOut().split(X, y)
        elif args.all:
            val_folds = [(np.array(range(len(X))), np.array([], dtype=int))]
            #make all/nothing split and also fix argparser
        
        y_pred = np.zeros(len(y)) - 1
        y_true = np.zeros(len(y)) - 1
        count=0
        for train_index, test_index in val_folds:
            print("fold: %s"%count)
            count += 1

            X_train, X_test = X[train_index], X[test_index]
            y_train, y_test = y[train_index], y[test_index]

            sampler = SMOTE(sampling_strategy='all')
            X_res, y_res = sampler.fit_resample(X_train, y_train)
            new_ind = np.random.permutation(range(len(y_res)))
            X_res = X_res[new_ind,:]
            y_res = y_res[new_ind]
#            if count < 2:
#                np.save(DESTDIR + "/yres", np.array(y_res))
            clf = RandomForestClassifier(n_estimators=n_est)
            clf.fit(X_res,y_res)
            if not args.all: # test model success
                y_pred[test_index] = clf.predict(X_test)
                y_true[test_index] = y_test
        if not args.all: #save model test results
            np.save(name_extension + '_ytrue', y_true)
            np.save(name_extension + '_ypred', y_pred)
            plot_confusion_matrix(y, y_pred, name_extension=name_extension + '_conf_matrix')
            importances = clf.feature_importances_
            importances = np.array(importances)
            np.save(DESTDIR + "/importances" + ext, importances, allow_pickle=True, fix_imports=True)
            plot_importances(importances, names, DESTDIR + "/importances%s.png" % ext)
        else: # save trained algorithm
            filename = DESTDIR + "/final_trained_rf%s.pkl" \
                    % ("_ia_only" if args.ia_only else "")
            pickle.dump(clf, open(filename, 'w+b'))
            
def main():
    start = time.time()
    
    if (args.kfold and args.all) or (args.all and args.loo)\
     or (args.loo and args.kfold) or not (args.all or args.kfold or args.loo):
         raise Exception("must use exactly one of the following flags: \
                         --kfold, --loo --all")

    X1, y1 = collect(1)  
    y1 = np.array(y1)
    
    run(X1, y1, NUM_TREES, DESTDIR + '/' + ext)

    end = time.time()
    print('running time: %s' % (end - start))

if __name__=='__main__':
    main()
    
    '''   
    
                                    Hi Hi Hi Hi Hi Hi
                             Hi.                               Hi
                          Hi.          Hi.       Hi.           Hi
                      Hi hi.                                       Hi hi
                      Hi hi.                Hi                    Hi. Hi
                           Hi.          HiHiHiHi.             Hi 
                                Hi.                               Hi 
                                      Hi hi hi hi hi hi hi 
                      Hi                       Hi
                      Hi                       Hi
                           Hi.                  Hi
                               Hi. Hi hi hi hi. Hi hi Hi
                                                 Hi.            Hi
                                                 Hi.          Hi
                                                 Hi.  Hi hi   
                                                 Hi
                                             Hi.   Hi 
                                          Hi.         Hi
                                       Hi.               Hi
    
    
    '''