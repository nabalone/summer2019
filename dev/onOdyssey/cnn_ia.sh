#!/bin/bash
#
#SBATCH --job-name=ia
#SBATCH --output=cnn_ia.log
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --time=600:00
#SBATCH --mem-per-cpu=20000

cd /n/holylfs03/LABS/berger_lab/nchou 
python cnn_aardvark_aug_concat.py -n 10 -a -b 50


###END OF THE FILE#####


