#!/bin/bash
#
#SBATCH --job-name=cnn_aardvark25weighted
#SBATCH --output=cnn_aardvark25weighted.log
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --time=1000:00
#SBATCH --mem-per-cpu=20000

cd /n/holylfs03/LABS/berger_lab/nchou 
python cnn_aardvark_weighted.py


###END OF THE FILE#####


