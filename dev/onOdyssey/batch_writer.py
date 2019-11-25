#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 23:28:49 2019

@author: nabalone
"""

def write_batch(flag_string):
    name = flag_string.replace('-', '_')
    name = name.replace(' ', '')
    string = '#!/bin/bash\n#\n#SBATCH --job-name=%s\n#SBATCH --output=cnn_%s.log\n#SBATCH -p shared\n#SBATCH --ntasks=1\n#SBATCH --time=600:00\n#SBATCH --mem-per-cpu=20000\n\ncd /n/holylfs03/LABS/berger_lab/nchou \npython cnn_aardvark_aug_concat.py %s \n\n###END OF THE FILE#####\n\n\n' %(name, name, flag_string)
    with open('batchfiles/cnn_run_%s.sh' % name, 'w+') as f:
        f.write(string)
    #!/bin/bash\n#\n#SBATCH --job-name=ia\n#SBATCH --output=cnn_ia.log\n#SBATCH -p shared\n#SBATCH --ntasks=1\n#SBATCH --time=600:00\n#SBATCH --mem-per-cpu=20000\n\ncd /n/holylfs03/LABS/berger_lab/nchou \npython cnn_aardvark_aug_concat.py -n 10 -a -b 50\n\n###END OF THE FILE#####\n\n\n'
    
for num in ['-a', '-3', '']:
    for conv in ['', '-c']:
        for batch in ['', '-b 200', 'b 50']:
            for arch in ['', '-p 100', '-d 0.5', '-d 0.25', '-p 10', '-p 10 -d 0.25']:
                flag_string = '%s %s %s %s' % (num, conv, batch, arch)
                write_batch(flag_string)