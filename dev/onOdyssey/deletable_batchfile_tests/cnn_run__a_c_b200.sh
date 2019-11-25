#!/bin/bash
#
#SBATCH --job-name=_a_c_b200
#SBATCH --output=cnn__a_c_b200.log
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --time=600:00
#SBATCH --mem-per-cpu=20000

cd /n/holylfs03/LABS/berger_lab/nchou 
python cnn_aardvark_aug_concat.py -a -c -b 200  

###END OF THE FILE#####


