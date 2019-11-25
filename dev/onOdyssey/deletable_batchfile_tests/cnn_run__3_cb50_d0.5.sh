#!/bin/bash
#
#SBATCH --job-name=_3_cb50_d0.5
#SBATCH --output=cnn__3_cb50_d0.5.log
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --time=600:00
#SBATCH --mem-per-cpu=20000

cd /n/holylfs03/LABS/berger_lab/nchou 
python cnn_aardvark_aug_concat.py -3 -c b 50 -d 0.5 

###END OF THE FILE#####


