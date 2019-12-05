#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 23:28:49 2019

@author: nabalone
"""
#TODO REMEMBER TO FIX MEM PER BATCH SIZE
def write_batch(flag_string):

    mem = 20000

#    if '200' in batch:
#        mem = 20000
#    elif '50' in batch:
#        mem = 15000
#    else:
#        mem = 30000

    name = flag_string.replace('-', '')
    name = name.replace(' ', '_')
    string = '#!/bin/bash\n#\n#SBATCH --job-name=%s\n#SBATCH --output=cnn_run_%s.log\n#SBATCH -p shared\n#SBATCH --ntasks=1\n#SBATCH --time=600:00\n#SBATCH --mem-per-cpu=%s\n\ncd /n/home00/nchou \npython cnn_aardvark_aug_concat.py %s \n\n###END OF THE FILE#####\n\n\n' %(name, name, mem, flag_string)
    with open('batchfiles/cnn_run_%s.sh' % name, 'w+') as f:
        f.write(string)
    #!/bin/bash\n#\n#SBATCH --job-name=ia\n#SBATCH --output=cnn_ia.log\n#SBATCH -p shared\n#SBATCH --ntasks=1\n#SBATCH --time=600:00\n#SBATCH --mem-per-cpu=20000\n\ncd /n/holylfs03/LABS/berger_lab/nchou \npython cnn_aardvark_aug_concat.py -n 10 -a -b 50\n\n###END OF THE FILE#####\n\n\n'


for batch in ['-b 50', '-b 199']:
    for mask in ['', '--mask']:
        for drop in ['-d 0.25', '']:
            for alt in ['', '--alt']:
                flag_string = '%s %s %s %s -c --mp 5 -n 300' % (batch, mask, drop, alt)
                write_batch(flag_string) 
#for num in ['']:#, '-3', '']:
#for lr in ['-l 0.00005', '-l 0.000005', '']:
#        for batch in ['', '-b 174', '-b 58']:
#            for arch in ['', '-p 5', '--mp 5', '-d 0.25', '-c', '-c --mp 5 -d 0.25', '-c --mp 5']:
#                    for mask in ['--mask', '']:
#                        flag_string = '%s %s %s %s -n 200' % (lr, batch, arch, mask)
#                        write_batch(flag_string, batch)
