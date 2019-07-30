import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.io import fits
from astropy import wcs
from astropy.nddata.utils import Cutout2D
from astropy.coordinates import SkyCoord
import astropy.units as u

#Read in my SNe file
#sne = np.loadtxt('./ps1confirmed_only_sne_without_outlier.txt',
#	usecols=0,dtype=str,skiprows=1)
#Read in Alerts table & look for SN
db = pd.read_table('alertstable_v3',sep=None,index_col = False)
db2 = pd.read_table('alertstable_v3.lasthalf',sep=None,index_col = False)
db = db.append(db2,ignore_index=True)
files = ['md06s072.g.stack_62.notyr1.fits']

#open up correct fits file
for sn in sne:
	sn = sn[3:]
	gdb = db.where(db['eventID'] == int(sn)).dropna()
	if len(gdb)>0:
		yr = gdb['year'].values[0]
		fname = str(gdb['field'].values[0])+'.g.stack_'+str(gdb['amp'].values[0])[:-2]+'.notyr1.fits'
		if fname in files:
			ra = gdb['ra'].values[0]
			dec = gdb['dec'].values[0]
			c = SkyCoord(ra+dec,unit=(u.hourangle,u.deg))

			f = fits.open(fname)
			w = wcs.WCS(f[0].header)
			newf = fits.PrimaryHDU()
			newf.data = f[0].data
			newf = Cutout2D(newf.data,c,size=[240,240],wcs=w)
			f[0].header['CRPIX1'] = newf.wcs.wcs.crpix[0]
			f[0].header['CRPIX2'] = newf.wcs.wcs.crpix[1]
			hdu = fits.PrimaryHDU(data=newf.data, header = f[0].header)
			hdu.writeto('new.fits')


#do a cut around the transient position of 1 arcmin


#save as new file