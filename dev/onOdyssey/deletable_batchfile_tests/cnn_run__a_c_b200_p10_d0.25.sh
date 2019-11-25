#!/bin/bash
#
#SBATCH --job-name=_a_c_b200_p10_d0.25
#SBATCH --output=cnn__a_c_b200_p10_d0.25.log
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --time=600:00
#SBATCH --mem-per-cpu=20000

cd /n/holylfs03/LABS/berger_lab/nchou 
python cnn_aardvark_aug_concat.py -a -c -b 200 -p 10 -d 0.25 

###END OF THE FILE#####


