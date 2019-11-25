#!/bin/bash
#
#SBATCH --job-name=_3_c_d0.25
#SBATCH --output=cnn__3_c_d0.25.log
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --time=600:00
#SBATCH --mem-per-cpu=20000

cd /n/holylfs03/LABS/berger_lab/nchou 
python cnn_aardvark_aug_concat.py -3 -c  -d 0.25 

###END OF THE FILE#####


