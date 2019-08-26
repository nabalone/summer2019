import glob
import os
SOURCEDIR = "C:/Users/Faith/Desktop/noey2019summer/ps1hosts"
DESTDIR = os.getcwd() + '/deblend0.05'

filenames = sorted(glob.glob(SOURCEDIR + '/psc*'))

import random
with open('presetradomfiles.txt', 'w+') as f:
    for i in range(100):
        f.write('"')
        f.write(filenames[int(len(filenames)*random.random())])
        f.write('", \n')