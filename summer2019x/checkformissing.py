# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 13:54:08 2019

@author: Faith
"""

import glob
SOURCEDIR = "F:\ps1hosts"
filenames = sorted(glob.glob(SOURCEDIR + '\\psc*.fits'))

for i in range(len(filenames)-1):
    if not ((filenames[i][22:23] == '6' and filenames[i+1][22:23] =='3')
            or (int(filenames[i+1][22:23]) == 
               int(filenames[i][22:23]) + 1)):

        print(filenames[i])
        print(filenames[i+1])
        print('')
        
    