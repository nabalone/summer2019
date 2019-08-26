#!/bin/bash
#
#SBATCH --job-name=cnn_aardvark25
#SBATCH --output=cnn_aardvark25.log
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --time=1000:00
#SBATCH --mem-per-cpu=20000

cd /n/holylfs03/LABS/berger_lab/nchou 
python cnn_aardvark.py


###END OF THE FILE#####


