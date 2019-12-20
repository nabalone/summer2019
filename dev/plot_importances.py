# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 12:12:50 2019

@author: Faith
"""
import numpy as np
import matplotlib.pyplot as plt
DIR = '/mnt/c/Users/Noel/Desktop/summer2019/dev/apple_run'
aa = np.load(DIR + "/importances3.npy")
a = []
for j in range(10):
    for i in range(4):
        if aa[10*i+j]== 0:
            continue
        a.append(aa[10*i+j])
#a.append(aa[40])
print(a)
colors = []
for i in ['r','b','g','c','m','k']:
    for j in range(4):
        colors.append(i)
colors.extend(colors[:])
plt.bar(range(36), a, color=colors[:36])
plt.plot(np.arange(-1,38,1),[aa[40]]*39, linestyle = "--", color='k')

plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off

plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off

plt.locator_params(axis = 'y', nbins = 5)
plt.gcf().subplots_adjust(left=0.3)
plt.savefig(DIR + "/importances_good.png")
