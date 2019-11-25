#!/bin/bash
#
#SBATCH --job-name=_d0.5
#SBATCH --output=cnn__d0.5.log
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --time=600:00
#SBATCH --mem-per-cpu=20000

cd /n/holylfs03/LABS/berger_lab/nchou 
python cnn_aardvark_aug_concat.py    -d 0.5 

###END OF THE FILE#####


