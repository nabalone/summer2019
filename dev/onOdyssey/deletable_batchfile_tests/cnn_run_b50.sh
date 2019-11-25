#!/bin/bash
#
#SBATCH --job-name=b50
#SBATCH --output=cnn_b50.log
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --time=600:00
#SBATCH --mem-per-cpu=20000

cd /n/holylfs03/LABS/berger_lab/nchou 
python cnn_aardvark_aug_concat.py   b 50  

###END OF THE FILE#####


