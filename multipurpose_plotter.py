# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 13:35:10 2019

@author: Noel
"""
import csv

X_PROP = 'Z'
Y_PROP = 'Abs_mag'

TYPES = ['SNIIn', 'SNIa', 'SNII', 'SNIbc', 'SLSNe']
TYPE_COLORS = {'SNIIn':'co', 'SNIa':'ro', 'SNII': 'bo', 'SNIbc':'go', 'SLSNe': 'mo'}

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

with open(CSVFILE, 'r') as f:
    r = csv.reader(f)
    r.readrow()
    r.readrow()
    for row in r:
        idNum = row[0]
        sn_type = typeDict[pad(idNum)]
        x_prop[sn_type].append(row[x_prop_num])
        y_prop[sn_type].append(row[y_prop_num])
        

areas = {}
# plot magnitude vs. area
for snType in TYPES:
    plt.plot(x_prop[snType], y_prop[snType], TYPE_COLORS[snType])
    
    plt.xlabel(X_PROP)
    plt.ylabel(Y_PROP)
    plt.savefig(PLOT_DIR + X_PROP + '_vs_' + Y_PROP, dpi=150)
    plt.show()
    plt.close()