#!/bin/bash
#
#SBATCH --job-name=_ab50_p100
#SBATCH --output=cnn__ab50_p100.log
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --time=600:00
#SBATCH --mem-per-cpu=20000

cd /n/holylfs03/LABS/berger_lab/nchou 
python cnn_aardvark_aug_concat.py -a  b 50 -p 100 

###END OF THE FILE#####


