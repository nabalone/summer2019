#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 21 16:46:28 2020

@author: nabalone
"""

import numpy as np
import matplotlib.pyplot as plt
import glob
import os

#PROJ_HOME = os.environ['DATA_SRCDIR'] #base of repo
#OUTPUTDIR = SOURCEDIR + '/outputs'
OUTPUTDIR = 'outputs'
predfiles = glob.glob(OUTPUTDIR + '/cnn_kfold_results/y_pred*.npy')
#trues = [pred.replace('pred', 'true') for pred in preds]
ypred_all = []
ytrue_all = []
for p in predfiles:
    ypred = np.load(p)
    ytrue = np.load(p.replace("pred", "true"))
    ypred_all.append(ypred)
    ytrue_all.append(ytrue)
    
ypred_all = np.hstack(ypred_all)
ytrue_all=np.hstack(ytrue_all)
print(ytrue_all.shape)
print(ypred_all.shape)



