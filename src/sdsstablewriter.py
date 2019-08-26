# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 16:55:35 2019

@author: Faith
"""
import numpy as np
import pandas as pd
import os
import glob
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
import os

start = time.time()
errs= []
SOURCEDIR = os.getcwd() + "/ps1hosts"

filenames = sorted(glob.glob(SOURCEDIR + '/psc*.[3-6].fits'))
fileset = set()
for f in filenames:
    dotSplit = f.split('.')
    idNumString = dotSplit[-3].split('c')[-1]
    fileset.add(idNumString)
    
print(len(fileset))
    
db = pd.read_table('alertstable_v3',sep=None,index_col = False, 
               engine='python')
db2 = pd.read_table('alertstable_v3.lasthalf',sep=None,index_col = False, 
                engine='python')
db = db.append(db2,ignore_index=True)
    
full_table = None    
index = 0
eventdict = {}
for f in sorted(list(fileset)):
    filename = glob.glob(SOURCEDIR + '/psc' + f + '.[3-6].fits')[0]
    image_file = fits.open(filename)
    
    # get event pixel coords
    w = WCS(filename)
    event = db.where(db['eventID'] == int(f)).dropna()
    eventRa = event['ra'].values[0] #values gives np arrays
    eventDec = event['dec'].values[0] #'hh:mm:ss.sss'
        # converting to degrees
    eventCoords = SkyCoord(eventRa, eventDec, unit=(u.hourangle, u.deg))
    eventRa = eventCoords.ra.deg
    eventDec = eventCoords.dec.deg
    
    eventX, eventY = w.all_world2pix(eventRa, eventDec, 1)
    

    # to get image dimensions in wcs:
    maxX = image_file[0].header['NAXIS1']
    maxY = image_file[0].header['NAXIS2']
    maxRa, maxDec = w.all_pix2world(1,maxY,1)
    minRa, minDec = w.all_pix2world(maxX,1,1)

    # fix formatting
    maxRa = maxRa.item(0)*u.deg
    minRa = minRa.item(0)*u.deg
    maxDec = maxDec.item(0)*u.deg
    minDec = minDec.item(0)*u.deg
        
    
    
    # make query
    print(f)
    query = "SELECT p.ra, p.dec, p.type, p.modelMag_g, p.modelMag_r, \
                p.modelMag_i, p.modelMag_z, pz.z, pz.zErr\
            FROM Photoz AS pz RIGHT JOIN PhotoObj AS p ON pz.objid = p.objid \
            WHERE p.mode = 1 AND p.ra < %s and p.ra > %s AND p.dec < %s and p.dec > %s" \
            % (maxRa, minRa, maxDec, minDec)
    sdssTable = SDSS.query_sql(query)
    if not sdssTable:
        continue
    sdssTable['idnum'] = [f]*len(sdssTable)
    eventdict[f] = index
    index += 1
    if not full_table:
        full_table = sdssTable
    else:    
        try:
            full_table = vstack([full_table, sdssTable])
        except Exception as e: 
            sdssTable.replace_column('z', [0]*len(sdssTable))
            sdssTable.replace_column('zErr', [-999.]*len(sdssTable))
            full_table = vstack([full_table, sdssTable])
            #raise
            #errs.append(str(f) + str(e))
    
#full_table = full_table.group_by('idnum')

full_table.write('sdss_queries.dat', format='ascii', overwrite=True)
'sdss_queries.dat'

with open('sdss_queries_index.txt', 'w+') as indexfile:
    indexfile.write(str(eventdict))
    
end = time.time()
print(end - start)