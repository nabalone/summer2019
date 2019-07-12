# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 13:35:10 2019

@author: Noel
"""
import csv
import numpy as np
import matplotlib.pyplot as plt
'''SLSNe:
330114
340195
440420
090022
060270
180279
120031
160103
480552
150381
110405
030129
590123
'''
import os

PLOT_DIR = os.getcwd() + '/msc_plots/'
if not os.path.isdir(PLOT_DIR):
    os.mkdir(PLOT_DIR)
    
CSVFILE = os.getcwd() + "/goodFourthRun/galaxiesdata_manuallycut.csv"

HEADER =['ID']
#perImageHeaders = ['KronRad', 'separation', 'x', 'y', 'RA', 'DEC', 'KronMag', 'Angle', 'Ellipticity']
perImageHeaders = ['KronRad (kpc)', 'separation (kpc)', 'area (kpc^2)', 'sep/area (kpc)',
                   'x', 'y','KronMag', 'Abs. Mag', 'Angle',
                   'Ellipticity', 'RA', 'Host RA', 'DEC', 'Host Dec',
                   'Discrepency (arcsecs)', 'Z', 'SDSS Photoz', 'pixelRank', 'chance coincidence']
for i in range(3,7):
    for val in perImageHeaders:
        HEADER.append(val + '_' + str(i))
HEADER = np.array(HEADER)

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

def run(X_PROP, Y_PROP):
    
    x_prop_num = np.where(HEADER==X_PROP)[0][0]
    y_prop_num = np.where(HEADER==Y_PROP)[0][0]

    with open(CSVFILE, 'r') as f:
        r = csv.reader(f)
        r.next()
        r.next()
        x_prop = {}
        y_prop = {}
        for t in TYPES:
            x_prop[t] = []
            y_prop[t] = []
        for row in r:
            idNum = row[0]
            sn_type = typeDict[pad(idNum)]
            try:
                float(row[x_prop_num])
                float(row[y_prop_num])
            except:
                continue
            x_prop[sn_type].append(float(row[x_prop_num]))
            y_prop[sn_type].append(float(row[y_prop_num]))
            
    def trim(string):
        fixed = string.replace(" ", "")
        fixed = fixed.replace("(", "")
        fixed = fixed.replace(")", "")
        fixed = fixed.replace("/", "Per")
        return fixed
    
    areas = {}
    # plot magnitude vs. area
    plot_args = []
    for snType in TYPES:
        plt.plot(x_prop[snType], y_prop[snType], TYPE_COLORS[snType], label=snType)
    #plt.plot(*plot_args)
    plt.xlabel(X_PROP)
    plt.ylabel(Y_PROP)
    plt.yscale("log")
    plt.legend()
    plt.savefig(PLOT_DIR + trim(X_PROP) + '_vs_' + trim(Y_PROP) + ".png", dpi=150)
    plt.show()
    plt.close()
    
run('KronRad (kpc)_3', 'separation (kpc)_3')
run('area (kpc^2)_4', 'KronMag_4')
run('sep/area (kpc)_5',  'Abs. Mag_5')
run('Ellipticity_6', 'Z_6')
run('pixelRank_4',  'chance coincidence_4')
run('SDSS Photoz_5', 'Discrepency (arcsecs)_5')