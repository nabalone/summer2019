# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 15:02:48 2019

@author: Faith
"""
import numpy as np
import pandas as pd
import os
import glob
import random
print(1)
#import csv
#import sep
#import ast
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
#from matplotlib.patches import Ellipse, RegularPolygon
from astropy.io import fits
from astropy.wcs import WCS
print(2)
from astropy import units as u
from astropy.coordinates import SkyCoord
#from astroquery.sdss import SDSS
#from astropy.cosmology import Planck13 as cosmo
#from astropy.table import Table, vstack
#from sdssTableGroupIndex import sdssTableGroupIndex

cols1 = ['KronRad (kpc)_3', 'separation (kpc)_3', 'area (kpc^2)_3', 'sep/sqrt(area) (kpc)_3', \
 'KronMag_3', 'Abs. Mag_3','Ellipticity_3', 'Z_3', 'pixelRank_3', 'chance coincidence_3']
cols = ['ID']
cols.extend(cols1[:])
for i in range(4,7):
    for e in cols1:
        cols.append(e[:-1] + str(i))

print(3)
test_indices = [37,
 38,
 39,
 40,
 153,
 154,
 155,
 156,
 268,
 269,
 270,
 286,
 287,
 288,
 305,
 306,
 307,
 308,
 366,
 367,
 368]

SOURCEDIR = os.getcwd() 
BASEDIR = "../../src/"
MASKDIR = BASEDIR
PIXDIR = BASEDIR + "combined/"#"ps1hosts/"#../vvillar/ps1_host/fits/"
CSVFILE = '../old_results/goodSeventhRun/galaxiesdata2.csv' # os.getcwd() + '/goodSeventhRun/untouched.csv'

filenames = sorted(glob.glob(PIXDIR + 'psc*.[3].fits'))
x_train = []
y_train = []
x_test = []
y_test = []

'''load event type dictionary'''
typeDict = {}
typefile = open(BASEDIR + 'ps1confirmed_only_sne_without_outlier.txt', 'r')
typefile.readline() #get rid of header
for line in typefile:
    parts = line.split()
    eventId = parts[0][3:]
    eventType = parts[1]
    typeDict[eventId] = eventType
typefile.close()
print(4)

zdict = {}
zfile = open(BASEDIR + 'new_ps1z.dat', 'r')
zfile.readline() #get rid of header

for line in zfile:
    parts = line.split()
    eventId = parts[0][3:]
    redshift = float(parts[1])
    zdict[eventId] = redshift
zfile.close()
zdict['100014'] = 0.357
zdict['300220'] = 0.094
zdict['380108'] = 0.159

intDict= {'SNIa':0, 'SNIbc':1, 'SNII':2, 'SNIIn':3, 'SLSNe':4}

'''load event location dictionary'''
db = pd.read_table(BASEDIR + 'alertstable_v3',sep=None,index_col = False,
               engine='python')
db2 = pd.read_table(BASEDIR + 'alertstable_v3.lasthalf',sep=None,index_col = False,
                engine='python')
db = db.append(db2,ignore_index=True)

def size(axis, image_data):
    return int(np.shape(image_data)[axis])

def pad(axis, side, image_data):
    width = 240 - size(axis, image_data)
    if axis == 0:
        zeros = np.zeros((width, size(1, image_data)))
        if side == 'below':
            return np.vstack((zeros, image_data))
        elif side == 'above':
            return np.vstack((image_data, zeros))
    else:
        zeros = np.zeros((size(0, image_data), width))
        if side == 'below':
            return np.hstack((zeros, image_data))
        elif side == 'above':
            return np.hstack((image_data, zeros))
        

def plot(data, myEventX, myEventY, myVmin=None, myVmax=None):
    fig, ax = plt.subplots()

    if myVmax == None or myVmin == None:
        _im = ax.imshow(data, interpolation='nearest', cmap='gray')
    else:
        _im = ax.imshow(data, interpolation='nearest', cmap='gray',
                        vmin = myVmin, vmax = myVmax)
    # triangle on event location
    p = RegularPolygon((int(myEventX), int(myEventY)), 3, radius=3)
    p.set_edgecolor('purple')
    p.set_facecolor('purple')
    ax.add_artist(p)
    plt.show()
    plt.close()

#for filename in ['C:\Users\Faith\Desktop\NoeySummer2019\summer2019/ps1hosts\psc000137.3.fits']:#filenames[:1]:
def load_data():

    masks = np.load(MASKDIR + 'all_masks.npz')

    #x_train = []
    x_test = [[],[],[],[],[]]
    x_sep = [[],[],[],[],[]]
    data = pd.read_csv(CSVFILE)
    X = data.loc[:, cols].values
    X = np.nan_to_num(X)
    #y_train = []
    #y_test = []
    #randfile = open('randomly_ordered_filenames2.txt', 'r')
    filenames = glob.glob(PIXDIR + 'psc*.[3].fits')
    for full_filename in filenames:#index in test_indices:
        print(full_filename)
        #filename = filenames[index]
        #full_filename = full_filename[:-1] #drop ending \n
        dotSplit = full_filename.split('.')
        idNumString = dotSplit[-3].split('c')[-1] #used for getting hectospec data
        idNum  = int(idNumString)
        all_colors = []
        if len(glob.glob(full_filename[:-6] + '[3-6]' + full_filename[-5:])) !=4:
            print("missing filters for %s" % idNum)
            continue 
        for i in range(3,7):
            filename = full_filename[:-6] + str(i) + full_filename[-5:]   
            with fits.open(filename) as image_file:
                image_data = image_file[0].data
            
            # find pixel coordinates of event
            event = db.where(db['eventID'] == idNum).dropna()
            eventRa = event['ra'].values[0] #values gives np arrays
            eventDec = event['dec'].values[0] #'hh:mm:ss.sss'
            
            # converting to degrees
            eventCoords = SkyCoord(eventRa, eventDec, unit=(u.hourangle, u.deg))
            eventRa = eventCoords.ra.deg
            eventDec = eventCoords.dec.deg
            #TODO: ircs or fk5
            
            #plot(image_data, eventX, eventY)
            
            # get event pixel coords
            w = WCS(filename)
            eventX, eventY = w.all_world2pix(eventRa, eventDec, 1)
            val = image_data[int(eventX), int(eventY)]
            paddedX = None
            paddedY = None
            if size(0, image_data) != 240:
                if eventX < 120:
                    paddedX = 'above'
                    image_data = pad(0, 'above', image_data)
                else:
                    paddedX = 'below'
                    image_data = pad(0, 'below', image_data)
            if size(1, image_data) != 240:
                if eventY < 120:
                    paddedY = 'above'
                    image_data = pad(1, 'above', image_data)
                else:
                    paddedY = 'below'
                    image_data = pad(1, 'below', image_data)
                    
            if size(0, image_data) != 240 or size(1, image_data) != 240:
                print('width, height: %s, %s' % (size(0, image_data), size(1, image_data)))
                raise Exception("bad dimensions")
            all_colors.append(image_data)    
            #plot(image_data, eventX, eventY)
        shp = image_data.shape
        redshift = zdict[idNumString]
        all_colors.append(np.full(shp, redshift))

        #maskfiles = glob.glob('masks/mask_%s.npy'%idNumString)
        #if len(maskfiles) > 1:
        #    raise
        #mask = np.load(maskfiles[0])
        mask = masks[idNumString]
        if paddedX:
            mask = pad(0, paddedX, mask)
        if paddedY:
            mask = pad(1, paddedY, mask)
        if mask.shape != all_colors[0].shape:
            raise
        all_colors.append(mask)

        all_colors = np.nan_to_num(np.array(all_colors))
        typ = intDict[typeDict[idNumString]]
        x_test[typ].append(all_colors)

        #lookup IDnum in SEP properties table from csv
        #print(np.asarray(X==idNum))
        #print(X[0])
        #print(X[1][0])
        rows, columns = np.asarray(X==idNum).nonzero()
        added = False
        #print(rows)
        #print(columns)
        for i in range(len(rows)):
            if columns[i] == 0: # match is in the ID columnd
                props = X[rows[i]]
                x_sep[typ].append(props)
                added = True
                break
        if not added:
            print("no properties found for " + idNumString)

        #y_test.append(typ)

#        if len(y_test) < 200:
#            x_test.append(all_colors)
#            y_test.append(intDict[typeDict[idNumString]])
#        else:
#            x_train.append(all_colors)
#            y_train.append(intDict[typeDict[idNumString]])
    for i in range(len(x_test)):
        x_test[i] = np.array(x_test[i])
        x_test[i] = np.transpose(x_test[i], [0,2,3,1])
        np.save("x_all2_%s" % i, x_test[i])
        np.save("x_sep2_%s" % i, x_sep[i])
#    x_train = np.array(x_train)
#    x_train = np.transpose(x_train, [0,2,3,1])
    #y_test = np.array(y_test)
#    y_train = np.array(y_train)
    #randfile.close()
    #np.save("x_all", x_test)
    #np.save("y_all", y_test)
    #return(((x_train,y_train),(x_test,y_test)))

load_data()
