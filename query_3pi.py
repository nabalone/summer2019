from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from reproject import reproject_interp
import os
from astropy.wcs import WCS, _wcs
from astropy.visualization import PercentileInterval, ImageNormalize
import requests
from astropy.table import Table
from astropy.nddata import CCDData, NDData
from io import BytesIO
import scipy
import warnings

# Coordinates in Degrees
#RA = 149.7495973
#151.7120204
#DEC = 3.157574165
#1.6927538

def get_boundaries(RA, DEC):
    # Image size in arcmin
    search_diameter = 1
    search_diameter_pixels = int(search_diameter / 0.25 * 60)

    # Get box of edges
    ra_min = RA - search_diameter / 60
    ra_max = RA + search_diameter / 60
    dec_min = DEC - search_diameter / 60
    dec_max = DEC + search_diameter / 60

    # Output image shape (pixels)
    shape = (search_diameter_pixels, search_diameter_pixels)

    # Create a new WCS object.  The number of axes must be set
    # from the start
    w = WCS(naxis=2)

    # Vector properties may be set with Python lists, or Numpy arrays
    w.wcs.crpix = [shape[0] / 2, shape[1] / 2]
    w.wcs.cdelt = np.array([0.25 / 3600, 0.25 / 3600])
    w.wcs.crval = [RA, DEC]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.pc = [[-1,0],[0,1]]
    
    return (ra_min, ra_max, dec_min, dec_max, w, shape)

def download_ps1_image(filename, saveas=None):
    """
    Download image from PS1 and correct luptitudes back to a linear scale.
    Parameters
    ---------------
    filename : PS1 image filename (from `get_ps1_filename`)
    saveas   : Path to save template file (default: do not save)
    Output
    ---------------
    ccddata : CCDData format of data with WCS
    """
    res = requests.get('http://ps1images.stsci.edu' + filename)
    hdulist = fits.open(BytesIO(res.content))

    # Linearize from luptitudes
    boffset = hdulist[1].header['boffset']
    bsoften = hdulist[1].header['bsoften']
    data_linear = boffset + bsoften * 2 * np.sinh(hdulist[1].data * np.log(10.) / 2.5)
    warnings.simplefilter('ignore')  # ignore warnings from nonstandard PS1 header keywords
    ccddata = CCDData(data_linear, wcs=WCS(hdulist[1].header), unit='adu')

    # Save the template to file
    if saveas is not None:
        ccddata.write(saveas, overwrite=True)

    return ccddata

def download_references(ra_min, dec_min, ra_max, dec_max, mag_filter, template_basename=None, catalog=None):
    """
    Download 1 to 4 references from PS1 as necessary to cover full RA & dec range
    Parameters
    ---------------
    ra_min, ra_max   : Minimum and Maximum RA and DEC
    dec_min, dec_max   in units of degrees
    mag_filter       : Filter color 'g', 'r', 'i', 'z', or 'y'
    template_basename: Filename of the output(s), to be suffixed by 0.fits, 1.fits, ...
    catalog          : Catalog to which to align the reference image WCS (default: do not align)
    Output
    ---------------
    refdatas   : List of CCDData objects containing the reference images
    """

    filename0 = get_ps1_filename(ra_min, dec_min, mag_filter)
    filename1 = get_ps1_filename(ra_max, dec_max, mag_filter)
    filename2 = get_ps1_filename(ra_min, dec_max, mag_filter)
    filename3 = get_ps1_filename(ra_max, dec_min, mag_filter)

    filenames = {filename0, filename1, filename2, filename3}
    refdatas = []
    for i, fn in enumerate(filenames):
        if template_basename is not None:
            saveas = template_basename + '{:d}.fits'.format(i)
            print('downloading', saveas)
        else:
            saveas = None
            print('downloading', fn)
        refdata = download_ps1_image(fn, saveas)
        if catalog is not None:
            _, stars = make_psf(refdata, catalog)
            try:
                refine_wcs(refdata.wcs, stars, catalog)
            except _wcs.InvalidTransformError:
                print('WARNING: unable to refine wcs')
        refdatas.append(refdata)

    return refdatas

def get_ps1_filename(ra, dec, filt):
    """
    Download Image from PS1 and correct luptitudes back to a linear scale.
    Parameters
    ---------------
    ra, dec : Coordinates in degrees
    filt    : Filter color 'g', 'r', 'i', 'z', or 'y'
    Output
    ---------------
    filename : PS1 image filename
    """

    # Query a center RA and DEC from PS1 in a specified color
    res = requests.get('http://ps1images.stsci.edu/cgi-bin/ps1filenames.py',
                 params={'ra': ra, 'dec': dec, 'filters': filt})
    t = Table.read(res.text, format='ascii')

    return t['filename'][0]

def assemble_reference(refdatas, wcs, shape):
    """Reproject and stack the reference images to match the science image"""
    refdatas_reprojected = []
    refdata_foot = np.zeros(shape, float)
    for data in refdatas:
        reprojected, foot = reproject_interp((data.data, data.wcs), wcs, shape)
        refdatas_reprojected.append(reprojected)
        refdata_foot += foot

    refdata_reproj = np.nanmean(refdatas_reprojected, axis=0)
    refdata_reproj[np.isnan(refdata_reproj)] = 0.
    refdata = CCDData(refdata_reproj, wcs=wcs, mask=refdata_foot == 0., unit='adu')
    return refdata

def run(RA, DEC, filter, template_filename): 
    ra_min, ra_max, dec_min, dec_max, w, shape = getBoundaries(RA, DEC)
    refdatas = download_references(ra_min, dec_min, ra_max, dec_max, filter)
    refdata = assemble_reference(refdatas, w, shape)
    refdata.write(template_filename, overwrite=True)
