# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 12:12:50 2019

@author: Faith
"""
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
PROJ_HOME = os.environ['DATA_SRCDIR']
sys.path.append(PROJ_HOME)
from src.random_forest_classifier import cols

def plot_importances(importances, feature_names, dest_name):
    a = []
    labels = []
    for j in range(10):
        for i in range(4):
            if importances[10*i+j]== 0:
                continue
            a.append(importances[10*i+j])
            labels.append(feature_names[10*i+j])
    #a.append(importances[40])
    colors = []
    for i in ['r','b','g','c','m','k']:
        for j in range(4):
            colors.append(i)
    colors.extend(colors[:])
    plt.bar(range(len(a)), a, color=colors[:len(a)])
    
    #plot line at height of random number input which is importances[40]
    plt.plot(np.arange(-1,len(a),1),[importances[40]]*(len(a)+1), linestyle = "--", color='k')
    plt.xticks(range(len(a)), labels, rotation='vertical', size=8)
    
    
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        ) 
    
    #plt.set_xticklabels(feature_names)
    
    plt.tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False) # labels along the bottom edge are off
    
    plt.locator_params(axis = 'y', nbins = 5)
    plt.gcf().subplots_adjust(left=0.3)
    plt.savefig(dest_name)
