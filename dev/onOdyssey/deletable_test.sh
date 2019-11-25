#!/bin/bash
#
#SBATCH --job-name=z
#SBATCH --output=cnn_z.log
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --time=600:00
#SBATCH --mem-per-cpu=20000

cd /n/holylfs03/LABS/berger_lab/nchou 
python cnn_aardvark_aug_concat.py -z 

###END OF THE FILE#####


