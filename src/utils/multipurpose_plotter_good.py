# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 13:35:10 2019

@author: Noel

"""



import matplotlib.pyplot as plt
import pandas as pd
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
PROJ_HOME = os.environ['DATA_SRCDIR']
OUTPUT_DIR = PROJ_HOME + '/src/outputs'
PLOT_DIR = PROJ_HOME + '/src/outputs'
if not os.path.isdir(PLOT_DIR):
    os.mkdir(PLOT_DIR)
    
CSVFILE = OUTPUT_DIR + '/galaxiesdata.csv'

TYPES = ['SNIIn', 'SNIa', 'SNII', 'SNIbc', 'SLSNe']
TYPE_COLORS = {'SNIIn':'co', 'SNIa':'ro', 'SNII': 'bo', 'SNIbc':'go', 'SLSNe': 'mo'}

COLORS = {'SNIa':'#d73027', 'SNIbc':'#fc8d59', 'SLSNe':'k',#'#fee090', 
          'SNII':'#91bfdb', 'SNIIn':'#4575b4'}
MARKERS = {'SNIIn':'s', 'SNIa':'*', 'SNII': 'v', 'SNIbc':'^', 'SLSNe': 'o'}

'''load event type dictionary'''
typeDict = {}
typefile = open(PROJ_HOME + '/src/ps1confirmed_added.txt', 'r')
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

#scatter plot Y_PROP vs X_PROP, which should both be column headers in the 
#galaxiesdata.csv file. Many of which end in underscore filter number.
# filter numbers 3-6 denote filters g,r,i,z in that order. 
# use dif_between=(filter1, filter2) to plot [Y_prop in filter1 - Y_prop in filter2] vs X_PROP
# e.g. to plot G-R color vs. redshift: run('redshift_3', 'Abs. Mag', dif_between=('g','r'))
#dif_between=(filternumber1, filternumber2) also works
#TODO why is redshift a filter-dependent property? this might be fixed, 
#in which case run('redshift', 'Abs. Mag', dif_between=('g','r'))
    
def run(X_PROP, Y_PROP, dif_between=None, log_scale_x=False, log_scale_y=False,
        ymax=None, ymin=None, xmax=None, xmin=None, xlabel=None, ylabel=None):#, plot_lim_mag=False):

    #r_prop_num = np.where(COLUMNS == 'Abs. Mag_4')[0][0]
    
    ptab = pd.read_table(OUTPUT_DIR + '/galaxiesdata.csv', sep=',', index_col=0)
    types = [typeDict[pad(i)] for i in ptab['ID']]
    ptab['type'] = types
    
    if X_PROP not in ptab.columns:
        raise Exception("Invalid X_PROP. Should be a column header in galaxiesdata.csv")        
    x_to_plot = {}
    for sn_type in TYPES:
        x_to_plot[sn_type] = ptab[ptab['type']==sn_type][X_PROP]
        
    y_to_plot = {}    
    if dif_between:
        #make sure we got 2 filters specified
        try:
            filt1, filt2 = dif_between
        except (TypeError, ValueError):
            raise Exception("dif_between argument should be None or a tuple of (filter1, filter2)")
            
        # remove filter specification from Y_PROP if different
        if Y_PROP[-2] == '_' and Y_PROP[-1] in ['3','4','5','6']:
            Y_PROP = Y_PROP[:-2]
            
        #convert to filter numbers if necessary
        filternum = {'g':3, 'r':4, 'i':5, 'z':6}
        if filt1 in ['g','r','i','z']:
            filt1 = filternum[filt1]
        if filt2 in ['g','r','i','z']:
            filt2 = filternum[filt2]
            
        if filt1 not in [3,4,5,6,'3','4','5','6'] or filt2 not in [3,4,5,6,'3','4','5','6']:
            raise Exception("Invalid filter specifications. dif_between should \
                            be a tuple of two filters, each specified by letter\
                            (g,r,i, or z) or number (3,4,5, or 6)")
            
        filt1_header = Y_PROP + '_' + str(filt1)
        filt2_header = Y_PROP + '_' + str(filt2)
        if filt1_header not in ptab.columns or filt2_header not in ptab.columns:
            raise Exception("Invalid Y_PROP. Should be a column header in galaxiesdata.csv")
        for sn_type in TYPES:
            y_to_plot[sn_type] = ptab[ptab['type']==sn_type][filt1_header] - ptab[ptab['type']==sn_type][filt2_header]
        
    else:
        if Y_PROP not in ptab.columns:
            raise Exception("Invalid Y_PROP. Should be a column header in galaxiesdata.csv")
        for sn_type in TYPES:
            y_to_plot[sn_type] = ptab[ptab['type']==sn_type][Y_PROP]
             
    def trim(string):
        fixed = string.replace(" ", "")
        fixed = fixed.replace("(", "")
        fixed = fixed.replace(")", "")
        fixed = fixed.replace("/", "Per")
        return fixed
    
#    if plot_lim_mag:
#        z= np.arange(0., 2., 0.1)
#        dL = cosmo.luminosity_distance(z)*1000/u.Mpc # in kpc
#        minMag = 25 - 5*np.log10(dL) - 10 + 2.5 * np.log10(1.+z)
#        plt.plot(z, minMag, label="limiting magnitude")
    
    for sn_type in TYPES:

        #y_prop[snType] = np.array(y_prop[snType])
        plt.plot(x_to_plot[sn_type], y_to_plot[sn_type], 
                 marker=MARKERS[sn_type], ms='5', linestyle="None", 
                 color=COLORS[sn_type], label=sn_type)
        #x_prop[snType]
    #plt.plot(*plot_args)
    
    font = {
        'weight' : 'normal',
        'size'   : 16}

    plt.rc('font', **font)
    plt.rc('figure', figsize=[10,7])
    
    
    if xlabel==None:
        xlabel=X_PROP
    if ylabel==None:
        if dif_between:
            ylabel = Y_PROP + " difference"
        else:
            ylabel = Y_PROP
            
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    
#TODO implement. Right now it doesn't work with log scale and axis limits, idk why
    if (log_scale_y and (ymin or ymax)) or (log_scale_x and (xmin or xmax)):
        raise Exception("log scale axis with axis min or max not yet implemented")
    
    if log_scale_x:
        plt.xscale("log")
    else:
        plt.xlim(xmax=xmax, xmin=xmin)
        
    if log_scale_y:
        plt.yscale("log")
    else:
        plt.ylim(ymax=ymax, ymin=ymin)
    #plt.legend(bbox_to_anchor=(1.1, 1.2))
    plt.legend()
    plt.tight_layout()
    plt.grid(b=True)
    #plt.gcf().subplots_adjust(bottom=0.15)
    plt.savefig(PLOT_DIR + '/' + trim(X_PROP) + '_vs_' + trim(Y_PROP) + ".png", dpi=150)
    plt.show()
    plt.close()
    
def main():  
    
    #run('Angle_6', 'separation (kpc)_3', log_scale_y=True)
    
    run('area (kpc^2)_4', 'KronMag_4', log_scale_x=True) 
    #run('redshift_3', 'Abs. Mag', dif_between=('g','r'), xlabel='Redshift', ylabel='G-R color', xmax=1, ymax=5)
    pass

#run('KronRad (kpc)_3', 'separation (kpc)_3')
#run('area (kpc^2)_4', 'KronMag_4')
#run('sep/area (kpc)_5',  'Abs. Mag_5')
#run('Ellipticity_6', 'Z_6')
#run('pixelRank_4',  'chance coincidence_4')
#run('SDSS Photoz_5', 'Discrepency (arcsecs)_5')

if __name__ == "__main__":
    main()