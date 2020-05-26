#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 10:40:52 2019

@author: nabalone
"""

import numpy as np
import matplotlib.pyplot as plt
import glob
import os

#TODO restore
PROJ_HOME = os.environ['DATA_SRCDIR'] #base of repo
SOURCEDIR = PROJ_HOME + '/src'
OUTPUTDIR = SOURCEDIR + '/outputs'
#OUTPUTDIR = 'outputs'


def load():
    predfiles = glob.glob(OUTPUTDIR + '/cnn_kfold_results/y_pred*.npy')
    ypred_all = []
    ytrue_all = []
    for p in predfiles:
        ypred = np.load(p)
        ytrue = np.load(p.replace("pred", "true"))
        ypred_all.append(ypred)
        ytrue_all.append(ytrue)
    global ypred_all
    global ytrue_all
    ypred_all = np.hstack(ypred_all)
    ytrue_all=np.hstack(ytrue_all)
    print(ytrue_all.shape)
    print(ypred_all.shape)

def plot_fom():
    pred = ypred_all
    test = ytrue_all
    actual_ias = []
    actual_ccs = []
    
    for i in range(len(test)):
        if test[i] == 0:
            actual_ias.append(pred[i][0])
        else:
            actual_ccs.append(pred[i][0])
    ntotia = len(actual_ias)
    fomx = []
    fomy = []
    effs = []
    pps = []
    for i in range(200):
        cutoff = i/100
        
        nsubia = 0
        for n in actual_ias:
            if n > cutoff:
                nsubia += 1
                
        nsubnonia = 0
        for n in actual_ccs:
            if n > cutoff:
                nsubnonia += 1
        if nsubia==0 and nsubnonia==0:
            continue
        eff = nsubia/ntotia
        pp = nsubia/(nsubia + 5*nsubnonia)
        effs.append(eff)
        pps.append(pp)
        fomy.append(eff*pp)
        fomx.append(cutoff)
        
    f = open(PROJ_HOME + '/src/foleymandel_fom_plus.txt')
    ff = f.readlines()
    xs = []
    ys = []
    exs = []
    eys = []
    pxs = []
    pys = []
    
    for i in ff:
        try:
            px, py, ex, ey, x, y = i.split()
        except:
            break
        xs.append(float(x))
        ys.append(float(y))
        exs.append(float(ex))
        eys.append(float(ey))
        pxs.append(float(px))
        pys.append(float(py))
    f.close()
    plt.plot(fomx, fomy, color='k', lw=3, label='Chou+')
    plt.plot(fomx, effs, color='r', lw=3, label='Chou+ efficiency')
    plt.plot(fomx, pps, color='m', lw=3, label='Chou+ psuedopurity')
    plt.plot(exs, eys, color='b', lw=3, label='Foley efficiency')
    plt.plot(pxs, pys, color='c', lw=3, label='Foley psuedopurity')
    plt.plot(xs, ys, color='g', lw=3, label='Foley-Mandel')
    plt.plot(np.arange(0,1.1,0.1),[np.max(ys)]*11, linestyle = "--", color='g')
    plt.plot(np.arange(0,1.1,0.1),[np.max(fomy)]*11, linestyle = "--", color='k')
    plt.axis([0, 1, 0, 1])
    plt.ylabel('Figure of Merit')
    plt.xlabel('Classification Weighting')
    plt.legend(loc=2)
    plt.savefig('fom_good_with_effandpp.png')
            
    