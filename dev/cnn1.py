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
import csv
import sep
import ast
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
from matplotlib.patches import Ellipse, RegularPolygon
from astropy.io import fits
from astropy.wcs import WCS

from astropy import units as u
from astropy.coordinates import SkyCoord
from astroquery.sdss import SDSS
from astropy.cosmology import Planck13 as cosmo
from astropy.table import Table, vstack
from sdssTableGroupIndex import sdssTableGroupIndex

filenames = sorted(glob.glob(PIXDIR + 'psc*.[3-6].fits'))
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
PIXDIR = SOURCEDIR + '/ps1hosts/'#"../vvillar/ps1_host/fits/"

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

'''load event location dictionary'''
db = pd.read_table('alertstable_v3',sep=None,index_col = False,
               engine='python')
db2 = pd.read_table('alertstable_v3.lasthalf',sep=None,index_col = False,
                engine='python')
db = db.append(db2,ignore_index=True)

def size(axis, image_data):
    return int(np.shape(image_data)[axis])

def pad(axis, side, image_data):
    width = 240 - size(axis, image_data)
    print width
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
for index in test_indices:
    filename = filenames[index]
    dotSplit = filename.split('.')
    idNumString = dotSplit[-3].split('c')[-1] #used for getting hectospec data
    idNum  = int(idNumString)
        
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
    print(val)
    if size(0, image_data) != 240:
        if eventX < 120:
            image_data = pad(0, 'above', image_data)
        else:
            image_data = pad(0, 'below', image_data)
    if size(1, image_data) != 240:
        if eventY < 120:
            image_data = pad(1, 'above', image_data)
        else:
            image_data = pad(1, 'below', image_data)
            
    if size(0, image_data) != 240 or size(1, image_data) != 240:
        print('width, height: %s, %s' % (size(0, image_data), size(1, image_data)))
        raise Exception("bad dimensions")
    
    #plot(image_data, eventX, eventY)