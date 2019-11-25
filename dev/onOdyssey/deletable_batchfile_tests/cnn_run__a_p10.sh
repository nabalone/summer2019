#!/bin/bash
#
#SBATCH --job-name=_a_p10
#SBATCH --output=cnn__a_p10.log
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --time=600:00
#SBATCH --mem-per-cpu=20000

cd /n/holylfs03/LABS/berger_lab/nchou 
python cnn_aardvark_aug_concat.py -a   -p 10 

###END OF THE FILE#####


