# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 16:55:35 2019

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
import sys
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
import time
import json

start = time.time()
errs= []
SOURCEDIR = "C:/Users/Faith/Desktop/noey2019summer/ps1hosts"

filenames = ['sn1bc_clean.csv', 'sniin_clean.csv', 'slsne_clean.csv']

idnumdict = json.load('idnumdict_3py.txt')
    
full_table = None    
index = 0
eventdict = {}
for filename in filenames:
    data = pd.read_csv(filename)
    type = file[:-4]
    for i in range(len(data)):
        eventRa = data['R.A.'][i].split(',')[0]
        eventDec = data['Dec.'][i].split(',')[0]
        z = data['z'][i]
        name = data['Name'][i]
    # converting to degrees
    eventCoords = SkyCoord(eventRa, eventDec, unit=(u.hourangle, u.deg))
    eventRa = eventCoords.ra.deg
    eventDec = eventCoords.dec.deg


    # fix formatting
    image_halfwidth = 0.00833 # = 0.5 arcmin
    maxRa = eventRa + image_halfwidth
    minRa = eventRa - image_halfwidth
    maxDec = eventDec + image_halfwidth
    minDec = eventDec - image_halfwidth
 
    # make query
    query = "SELECT p.ra, p.dec, p.type, p.modelMag_g, p.modelMag_r, \
                p.modelMag_i, p.modelMag_z, pz.z\
            FROM Photoz AS pz RIGHT JOIN PhotoObj AS p ON pz.objid = p.objid \
            WHERE p.mode = 1 AND p.ra < %s and p.ra > %s AND p.dec < %s and p.dec > %s" \
            % (maxRa, minRa, maxDec, minDec)
    sdssTable = SDSS.query_sql(query)
    if not sdssTable:
        continue
    sdssTable['idnum'] = [f]*len(sdssTable)
    idNum = idnumdict[name]
    eventdict[idNum] = index
    index += 1
    if not full_table:
        full_table = sdssTable
    else:    
        try:
            full_table = vstack([full_table, sdssTable])
        except Exception as e:
            sdssTable.replace_column('z', [0]*len(sdssTable))
            full_table = vstack([full_table, sdssTable])
            #raise
            #errs.append(str(f) + str(e))
    
#full_table = full_table.group_by('idnum')

full_table.write('sdss_queries.dat', format='ascii', overwrite=True)
'sdss_queries.dat'

with open('sdss_queries_index_psTo3pi.txt', 'w+') as indexfile:
    json.dump(eventdict, indexfile) #indexfile.write(str(eventdict))
    
end = time.time()
print end - start