#!/bin/bash
#
#SBATCH --job-name=_cb50_p10
#SBATCH --output=cnn__cb50_p10.log
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --time=600:00
#SBATCH --mem-per-cpu=20000

cd /n/holylfs03/LABS/berger_lab/nchou 
python cnn_aardvark_aug_concat.py  -c b 50 -p 10 

###END OF THE FILE#####


