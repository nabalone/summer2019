import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.io import fits
from astropy import wcs
from astropy.nddata.utils import Cutout2D
from astropy.coordinates import SkyCoord
import astropy.units as u
import glob

#Read in my SNe file
sne = np.loadtxt('./ps1confirmed_only_sne_without_outlier.txt',
    usecols=0,dtype=str,skiprows=1)
#Read in Alerts table & look for SN
db = pd.read_table('alertstable_v3',sep=None,index_col = False)
db2 = pd.read_table('alertstable_v3.lasthalf',sep=None,index_col = False)
db = db.append(db2,ignore_index=True)
#files = ['md06s072.g.stack_62.notyr1.fits']
jdir = '/n/holylfs03/LABS/berger_lab/djones01/missingtemplates/'
#jdir = '/n/holylfs03/LABS/berger_lab/djones01/v3tmpl/'
#open up correct fits file

sne = []
with open('wrong_year_still_missing.txt', 'r') as f:
    for line in f.readlines():
        sne.append(line.split()[0])
print(len(sne))
for sn in sne:
    colors = ['g','r','i','z']
    #sn = sn[3:]
    gdb = db.where(db['eventID'] == int(sn)).dropna()   
    for j,filt in enumerate(['3','4','5','6']):
        if len(gdb)>0:
            yr = str(int(gdb['year'].values[0]-8))
            if yr == '6':
                yr = '5'
            
            #fname = 'PSc'+sn+'.'+str(gdb['field'].values[0])+'.'+colors[j]+'.stack_'+str(gdb['amp'].values[0])[:-2]+'.notyr'+yr+'.sw.fits'
            
            fname = 'PSc'+sn+'.'+str(gdb['field'].values[0])+'.'+colors[j]+'.stack_'+str(gdb['amp'].values[0])[:-2]+'.notyr'+'*'+'.sw.fits'
            f1 = jdir+'PSc'+sn+'/tmpl/'+colors[j]+'/'
            fname=f1+fname
            fname = glob.glob(fname)[0]
            ra = gdb['ra'].values[0]
            dec = gdb['dec'].values[0]
            c = SkyCoord(ra+dec,unit=(u.hourangle,u.deg))
            try:
                f = fits.open(fname)
                print(fname)
            except Exception(e):
                print(e)
                print(sn,filt)
                continue
            w = wcs.WCS(f[0].header)
            newf = fits.PrimaryHDU()
            newf.data = f[0].data
            newf = Cutout2D(newf.data,c,size=[240,240],wcs=w)
            f[0].header['CRPIX1'] = newf.wcs.wcs.crpix[0]
            f[0].header['CRPIX2'] = newf.wcs.wcs.crpix[1]
            hdu = fits.PrimaryHDU(data=newf.data, header = f[0].header)
            hdu.writeto('from_any_year/'+ 'psc'+sn+'.'+filt+'.fits')
            with open('make_stamps_log.txt', 'a+') as logfile:
                logfile.write('%s %s' % (sn, filt))
            continue






            #NOISE
            fname = str(gdb['field'].values[0])+'.'+colors[j]+'.stack_'+str(gdb['amp'].values[0])[:-2]+'.notyr'+yr+'.noise.fits'
            f1 = jdir+'notyr'+yr+'/0x161'+filt+'/'
            fname=f1+fname
            ra = gdb['ra'].values[0]
            dec = gdb['dec'].values[0]
            c = SkyCoord(ra+dec,unit=(u.hourangle,u.deg))
            try:
                f = fits.open(fname)
            except:
                print(sn,filt)
                continue
            w = wcs.WCS(f[0].header)
            newf = fits.PrimaryHDU()
            newf.data = f[0].data
            newf = Cutout2D(newf.data,c,size=[240,240],wcs=w)
            f[0].header['CRPIX1'] = newf.wcs.wcs.crpix[0]
            f[0].header['CRPIX2'] = newf.wcs.wcs.crpix[1]
            hdu = fits.PrimaryHDU(data=newf.data, header = f[0].header)
            hdu.writeto('psc'+sn+'.'+filt+'.noise.fits')
            #MASK
            fname = str(gdb['field'].values[0])+'.'+colors[j]+'.stack_'+str(gdb['amp'].values[0])[:-2]+'.notyr'+yr+'.mask.fits'
            f1 = jdir+'notyr'+yr+'/0x161'+filt+'/'
            fname=f1+fname
            ra = gdb['ra'].values[0]
            dec = gdb['dec'].values[0]
            c = SkyCoord(ra+dec,unit=(u.hourangle,u.deg))
            try:
                f = fits.open(fname)
            except:
                print(sn,filt)
                continue
            w = wcs.WCS(f[0].header)
            newf = fits.PrimaryHDU()
            newf.data = f[0].data
            newf = Cutout2D(newf.data,c,size=[240,240],wcs=w)
            f[0].header['CRPIX1'] = newf.wcs.wcs.crpix[0]
            f[0].header['CRPIX2'] = newf.wcs.wcs.crpix[1]
            hdu = fits.PrimaryHDU(data=newf.data, header = f[0].header)
            hdu.writeto('psc'+sn+'.'+filt+'.mask.fits')
#do a cut around the transient position of 1 arcmin


#save as new file
