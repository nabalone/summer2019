# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 12:12:50 2019

@author: Faith
"""
import numpy as np
import matplotlib.pyplot as plt
a = np.load("sorted_order.npy")
colors = []
for i in ['r','b','g','c','m','k']:
    for j in range(4):
        colors.append(i)
colors.extend(colors[:])
plt.bar(range(42), a, color=colors[:42])
plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off

plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off

plt.locator_params(axis = 'y', nbins = 5)
plt.gcf().subplots_adjust(left=0.3)
plt.savefig("importances.png")