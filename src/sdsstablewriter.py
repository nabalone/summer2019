# -*- coding: utf-8 -*-
"""
Query SDSS database for all supernovae we have pictures for
to get photoz and types of nearby objects
Types are for disqualifying stars from being chosen hosts

Saves queries in 'sdss_queries.dat', indexed according to
'sdss_queries_index.txt'
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
import json
PROJ_HOME = os.environ['DATA_SRCDIR'] #base of repo

start = time.time()
SOURCEDIR = PROJ_HOME + '/src'
PIXDIR = SOURCEDIR + '/all_fits' #fits stamps
OUTPUTDIR = SOURCEDIR + '/outputs'

filenames = sorted(glob.glob(PIXDIR + '/psc*.[3-6].fits')) #list of fits files

#make a set of all the sn which have fits files in PIXDIR
sn_set = set()
for sn in filenames:
    dotSplit = sn.split('.')
    idNumString = dotSplit[-3].split('c')[-1] #6 digit 0-left-padded SN num from filename
    sn_set.add(idNumString)
    
print(len(sn_set))
    
db = pd.read_csv(SOURCEDIR + '/alertstable_v3',sep=None,index_col = False, 
               engine='python')
db2 = pd.read_csv(SOURCEDIR + '/alertstable_v3.lasthalf',sep=None,index_col = False, 
                engine='python')
db = db.append(db2,ignore_index=True)
    
full_table = None    
for sn in sorted(list(sn_set)):
    filename = glob.glob(PIXDIR + '/psc' + sn + '.[3-6].fits')[0]
    image_file = fits.open(filename)
    
    # find the ra and dec ranges covered by the image
    w = WCS(filename)
    maxX = image_file[0].header['NAXIS1'] #image dimentions in pixels
    maxY = image_file[0].header['NAXIS2']
    maxRa, maxDec = w.all_pix2world(1,maxY,1) #coordinates of the image corners
    minRa, minDec = w.all_pix2world(maxX,1,1)
        
    # make query
    query = "SELECT p.ra, p.dec, p.type, p.modelMag_g, p.modelMag_r, \
                p.modelMag_i, p.modelMag_z, pz.z, pz.zErr\
            FROM Photoz AS pz RIGHT JOIN PhotoObj AS p ON pz.objid = p.objid \
            WHERE p.mode = 1 AND p.ra < %s and p.ra > %s AND p.dec < %s and p.dec > %s" \
            % (maxRa, minRa, maxDec, minDec)
    #print(query)
    sdssTable = SDSS.query_sql(query) #query result
    if not sdssTable: #probably no objects in this image on sdss
        continue
    #record the corresponding supernova number in the query output table
    sdssTable['idnum'] = [sn]*len(sdssTable) 

    #append this query result to bottom of table
    if not full_table:
        full_table = sdssTable
    else:    
        try:
            full_table = vstack([full_table, sdssTable])
        except Exception as e:  #probably this object has not photoz in sdss
            #photozs are only used in host selection, with outliers expected
            sdssTable.replace_column('z', [-1]*len(sdssTable)) #so this should be okay
            sdssTable.replace_column('zErr', [-999.]*len(sdssTable))
            full_table = vstack([full_table, sdssTable])

#save the table of all the query results
full_table.write(OUTPUTDIR + '/sdss_queries.dat', format='ascii', overwrite=True)
    
end = time.time()
print(end - start)
